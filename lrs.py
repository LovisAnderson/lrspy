from gmpy2 import mpz, divexact
from lrs_datastructures import LrsDict, Variable
from copy import deepcopy
from abc import ABC, abstractmethod
from enum import Enum


class Lrs(ABC):
    def __init__(self, hyperplane_matrix, m, d, bounding_box=None):
        self.hyperplanes = deepcopy(hyperplane_matrix)
        self.matrix = hyperplane_matrix
        self.nr_hyperplanes = m
        self.B = LrsDict() # Basis
        self.C = LrsDict() # Cobasis
        self.m = m # Number of input hyperplanes
        self.d = d # dimension of embedding space + 1 (projective dimension)
        self.det = mpz(1) # determinant of the matrix, quasi the shared denominator
        self.cobases = [] # list of cobases found
        self.vertices = [] # list of vertices found
        self.position_vectors = [] # relative position to hyperplanes for each vertex
        self.i = self.d # Basis index in which we pivot
        self.j = 0 # Cobasis index in which we pivot
        self.boxed = True if bounding_box is not None else False # Flag that indicates if we pivot inside a box (given through constraints)
        self.bounding_box = bounding_box if bounding_box is not None else []
        self.drop_objective_value = True # Flag to indicate that we do not compute or update the target value of the objective (A[0}[d])

    def set_objective(self):
        for i in range(1, self.d):
            self.matrix[0][i] = mpz(-self.det)

    def augment_matrix_with_objective(self):
        objectiveRow = [mpz(1)]*(self.d)
        objectiveRow[0] = mpz(0)
        self.matrix.insert(0, objectiveRow)

    def init_basis(self):
        self.B = LrsDict()
        f = Variable(0)
        f.basis_index = 0
        f.box_variable = False
        f.slack_variable = False
        self.B.append(f)
        for i in range(self.d, self.d + self.m):
            b = Variable(i)
            b.box_variable = False
            b.slack_variable = True
            b.hyperplane_index = i - self.d
            self.B.append(b)
        self.B.order = list(range(self.m + 1))

    def init_cobasis(self):
        for i in range(1, self.d):
            c = Variable(i)
            c.box_variable = False
            c.slack_variable = False
            self.C.append(c)
        g = Variable(self.d + self.m)
        g.box_variable = False
        g.slack_variable = False
        self.C.append(g)
        self.C.order = list(range(1, self.d)) + [0]

    def init_dicts(self):
        self.init_basis()
        self.init_cobasis()

    def variables_into_basis(self):
        for k in range(self.d - 1):
            self.j = 0
            self.i = 1
            while (self.B[self.i] in range(1, self.d) or
                   self.matrix[self.B.order[self.i]][self.C.order[self.j]] == 0):
                self.i += 1
            self.pivot()
        print(self.info_string('Variables moved into basis.'))

    def first_basis(self):
        self.variables_into_basis()
        self.set_objective()
        if self.boxed:
            self.resort_inequalities() # todo necessary?
            self.move_into_box()
        self.make_feasible()
        self.resort_inequalities()
        self.append_solution()

    def resort_inequalities(self):
        # Sorts variables corresponding to inequalities s.t. they basis is 0, ..., m
        # Assumes initialized dicts with variables 0,..,d at start of basis and d+m at end of cobasis

        for i, b in enumerate(self.B[self.d:]):
            self.B[i + self.d] = self.B[i + self.d].change_variable(i + self.d)
        for i, c in enumerate(self.C[:-1]):
            self.C[i] = self.C[i].change_variable(self.m + 1 + i)

    def go_to_root(self):
        self.i, self.j = self.select_pivot()
        while self.i != 0 or self.j != 0:
            self.pivot()
            self.i, self.j = self.select_pivot()

    def make_feasible(self):
        for i, row in enumerate(self.matrix[self.d:self.m+1]):
            if row[0] < 0:
                self.matrix[i + self.d] = [-a for a in row]

    def move_into_box(self):

        while not self.inside_box():
            for i, b in enumerate(self.B):
                if not b.box_variable:
                    continue
                elif self.matrix[self.B.order[i]][0] < 0:
                    # Primal infeasible Variable
                    self.i = i
                    break
            while self.matrix[self.B.order[self.i]][self.C.order[self.j]] <= 0:
                # To get primal feasible variable the sign of A[Row[i]][0] has to change
                # This happens if A[Row[i]][Column[j]] > 0
                self.j += 1
            if self.j == self.d:
                raise ValueError
            self.pivot()
        print(self.info_string('After first basis with bounding box'))

    def add_box_constraints(self, constraints):
        """
        Dicts have to be initialized before.
        Constraints are of the form a0 + a1x1 + ... + a_d-1 x_d-1 >= 0
        """
        self.box_constraints = constraints
        self.matrix += constraints
        self.startBox = self.m + self.d
        box_variables = []
        for i in range(self.m + self.d, self.m + self.d + len(constraints)):
            box_v = Variable(i)
            box_v.box_variable = True
            box_v.slack_variable = True
            box_v.hyperplane_index = i - self.d
            box_variables.append(box_v)
        self.B += box_variables
        self.C[-1] = self.C[-1].change_variable(self.m + self.d + len(constraints))
        self.B.order += list(range(self.m + 1, self.m + 1 + len(constraints)))
        self.m += len(constraints)
        self.boxed = True

    def pivot_stays_in_box(self, i, j):
        pivotRow = self.B.order[i]
        pivotColumn = self.C.order[j]
        pivotElement = self.matrix[pivotRow][pivotColumn]
        insideBox = True

        for k, b in enumerate(self.B):
            if k == i:
                if not self.C[j].box_variable: # If a not box variable is pivoted in we do not care about sign
                    continue
                elif self.matrix_entry_after_pivot(self.B.order[k], 0, pivotRow, pivotColumn, pivotElement) < 0:
                    insideBox = False
                    break
            elif b.box_variable:
                # Determinant sign is changed before matrix update if pivotelement < 0
                # Therefore we need this multiplier to get True output
                detMultiplier = 1 if pivotElement > 0 else -1
                if detMultiplier * self.matrix_entry_after_pivot(
                    self.B.order[k], 0, pivotRow, pivotColumn, pivotElement
            ) < 0:
                    insideBox = False
                    break
        return insideBox

    @abstractmethod
    def select_pivot(self):
        pass

    def inside_box(self):
        for k, v in enumerate(self.B):
            if v.box_variable and self.matrix[self.B.order[k]][0] < 0:
                return False
        return True

    def search(self):
        logger.info(self.info_string('Search start:'))
        self.i = self.d
        nextbasis = True
        backtrack = False
        while nextbasis:
            self.j = 0
            self.i = self.d
            while self.j < self.d or self.B[self.m] != self.m:
                if self.j == self.d - 1 and self.B[self.m] == self.m:
                    logging.info('All bases found!')
                    logging.info('vertices', self.vertices)
                    nextbasis = False
                    yield SearchStatus.DONE
                    break
                if backtrack:
                    logging.debug('Pivoting back!')
                    self.i, self.j = self.select_pivot()
                    self.pivot()
                    self.increment()
                    logger.debug('i: {}, j: {}'.format(self.i, self.j))
                    backtrack = False
                    yield SearchStatus.BACKTRACKED
                else:
                    while self.j < self.d - 1 and not self.reverse():
                        self.increment()
                        yield SearchStatus.INCREMENTED
                    if self.j == self.d - 1:
                        backtrack = True
                    else:
                        if self.lex_min():
                            self.append_solution()
                        yield SearchStatus.NEWBASIS

                        logger.debug('start tree search from new root')
                        break

    def append_solution(self):
        logger.debug('Append basis: {}'.format(self.B))
        logger.debug('Vertex: {}'.format(self.get_vertex()))
        self.cobases.append(deepcopy(self.C))
        self.vertices.append(self.get_vertex())
        self.position_vectors.append(self.get_position_vector())

    def get_vertex(self):
        vertex = tuple(self.matrix[self.B.order[k]][0] / self.det for k in range(1, self.d))
        return vertex

    def get_position_vector(self):
        position_vector = [0]*self.nr_hyperplanes
        for i, b in enumerate(self.B):
            if b.slack_variable and not b.box_variable:
                position_vector[b.hyperplane_index] = 1 if self.matrix[self.B.order[i]][0] < 0 else -1
        return position_vector

    def update(self):
        B_out = deepcopy(self.B[self.i])
        C_Out = deepcopy(self.C[self.j])
        self.B[self.i], self.C[self.j] = self.C[self.j], self.B[self.i]
        self.B.sort_respecting_order()
        self.C.sort_respecting_order()
        self.i = self.B.index(C_Out)
        self.j = self.C.index(B_out)

    def reverse(self):
        logger.debug('In reverse: i: {}, j:{}'.format(self.i, self.j))
        possibleReversePivot = self.necessary_condition_for_reverse()
        if not possibleReversePivot:
            logger.debug('Not a possible reverse pivot!')
            return False
        self.pivot()
        i_forward, j_forward = self.select_pivot()
        if i_forward == self.i and j_forward == self.j:
            logger.debug('valid reverse')
            return True
        else:
            logger.debug('Not valid reverse: pivoting back')
            self.pivot()
            return False

    @abstractmethod
    def necessary_condition_for_reverse(self):
        pass

    def lex_min(self):
        return True

    def info_string(self, infoString=''):
        str = infoString

        str += 'Basis: {} \n'.format(self.B)
        str += 'Cobasis: {}\n'.format(self.C)
        str += 'det: {}\n'.format(format_int(self.det, with_plus_sign=False))
        str += 'matrix:\n'
        str += self.matrix_with_variables_str()
        return str

    def matrix_entry_after_pivot(self, i, j, pivotRow, pivotColumn, pivotElement):
        if self.drop_objective_value and i == 0 and j == 0:
            return mpz(0)
        if i == pivotRow:
            if j == pivotColumn:
                return self.det
            if pivotElement > 0:
                return self.matrix[i][j] * mpz(-1)
            else:
                return self.matrix[i][j]
        if j == pivotColumn:
            if pivotElement < 0:
                return self.matrix[i][j] * mpz(-1)
            else:
                return self.matrix[i][j]
        nominator = self.matrix[i][j] * pivotElement - self.matrix[i][pivotColumn] * self.matrix[pivotRow][j]
        return divexact(nominator, self.det)

    def pivot(self):
        logger.debug('pivot: outIndex: {}; inIndex: {}'.format(self.i, self.j))
        logger.debug('outVariable: {}; inVariable: {}'.format(self.B[self.i], self.C[self.j]))

        row = self.B.order[self.i]
        column = self.C.order[self.j]
        pivotElement = deepcopy(self.matrix[row][column])
        self.det = self.det if pivotElement > 0 else -self.det

        self.matrix = [
            [self.matrix_entry_after_pivot(i, j, row, column, pivotElement) for j in range(self.d)]
            for i in range(self.m +1)
        ]
        self.det = pivotElement if pivotElement > 0 else - pivotElement
        self.update()
        logger.debug(self.info_string('After pivot:'))

    def increment(self):
        if self.i == self.m:
            self.j += 1
            self.i = self.d
        else:
            self.i += 1

    def pretty_print_matrix(self):
        str = '[\n'
        for row in self.matrix:
            str += '['
            for entry in row:
                if entry.numerator == 0:
                    str += ' 0'
                elif entry.denominator == 1:
                    str += ' {}'.format(entry.numerator)
                else:
                    str += ' {}/{},'.format(entry.numerator, entry.denominator)
            str += ']\n'
        str += ']'
        print(str)

    def matrix_with_variables_str(self):
        str = ''
        for i, b in enumerate(self.B):
            row = []
            for j, c in enumerate(self.C):
                entry = format_int(self.matrix[self.B.order[i]][self.C.order[j]], with_plus_sign=False)
                row.append('A[{}][{}]={}, '.format(b, c, entry))
            str += '{:<12}  {:<12}  {:<12}'.format(*row)
            str += '\n'
        return str

    def hyperplane_variables(self):
        slack_variables = []
        for v in self.B + self.C:
            if not v.slack_variable:
                continue
            else:
                slack_variables.append(v)
        slack_variables.sort(key=lambda x: x.hyperplane_index)
        return slack_variables

    @staticmethod
    def max_index_of_smaller_number(list, number):
        """ Assumes sorted input list"""
        for i, element in enumerate(list):
            if element >= number:
                return i - 1


def hyperplane_string(hyperplane):
    hyperplane_str = format_int(hyperplane[0], with_plus_sign=False)
    for i, xi in enumerate(hyperplane[1:]):
        hyperplane_str += format_int(xi) + 'x<sub>{}</sub>'.format(i)
    hyperplane_str += '= 0'
    return hyperplane_str


def format_int(integer, with_plus_sign=True):
    if with_plus_sign:
        if abs(integer) < 1e+6:
            return format(integer, '+')
        else:
            return '{:+.3E}'.format(int(integer))
    else:
        if abs(integer) < 1e+6:
            return format(integer)
        else:
            return '{:.3E}'.format(int(integer))

class SearchStatus(Enum):
    NONE = "no status"
    BACKTRACKED = "backtracked"
    DONE = "all vertices found"
    INCREMENTED = "incremented"
    NEWBASIS = "new basis found"
