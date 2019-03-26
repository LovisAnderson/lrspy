from gmpy2 import mpz, divexact
import pytest
from copy import deepcopy
from abc import ABC, abstractmethod


class Lrs(ABC):
    def __init__(self, inequality_matrix, m, d):
        self.matrix = inequality_matrix
        self.inequality_ordering = [] # list that
        self.B = [] # Basis
        self.C = [] # Cobasis
        self.Row = [] # B[i] is to be found in Row[i] of matrix
        self.Column = [] # C[i] is to be found in Column[i] of matrix
        self.m = m # Number of input hyperplanes
        self.d = d # dimension of embedding space + 1 (projective dimension)
        self.det = mpz(1) # determinant of the matrix, quasi the shared denominator
        self.bases = [] # list of bases found
        self.vertices = [] # list of vertices found
        self.position_vectors = [] # relative position to hyperplanes for each vertex
        self.i = self.d # Basis index in which we pivot
        self.j = 0 # Cobasis index in which we pivot
        self.boxed = False # Flag that indicates if we pivot inside a box (given through constraints)

    def setObjective(self):
        for i in range(1, self.d):
            self.matrix[0][i] = mpz(-self.det)

    def augmentWithObjective(self):
        objectiveRow = [mpz(1)]*(self.d)
        objectiveRow[0] = mpz(0)
        self.matrix.insert(0, objectiveRow)

    def initDicts(self):
        self.B = [0] + list(range(self.d, self.d + self.m))
        self.inequality_ordering = list(range(self.d, self.d + self.m))
        self.Row = list(range(self.m + 1))
        self.C = list(range(1, self.d)) + [self.d + self.m]
        self.Column = list(range(1, self.d)) + [0]

    def firstBasis(self):
        for k in range(self.d - 1):
            self.j = 0
            self.i = 1
            while (self.B[self.i] in range(1, self.d) or
                   self.matrix[self.Row[self.i]][self.Column[self.j]] == 0):
                self.i += 1
            self.pivot()
        self.printInfo('After first basis')
        self.resort_inequalities()
        if not self.boxed:
            self.appendSolution()

    def resort_inequalities(self):
        # Sorts variables corresponding to inequalities s.t. they basis is 0, ..., m
        # Assumes initialized dicts with variables 0,..,d at start of basis and d+m at end of cobasis
        for i, b in enumerate(self.B[self.d:]):
            self.inequality_ordering[i] = b
            self.B[i + self.d] = i + self.d
            if self.boxed and b in self.boxIndices:
                self.boxIndices.remove(b)
                self.boxIndices.add(i)
        for i, c in enumerate(self.C[:-1]):
            self.inequality_ordering[self.m - self.d + 1 + i] = c
            self.C[i] = self.m + 1 + i
            if self.boxed and b in self.boxIndices:
                self.boxIndices.remove(b)
                self.boxIndices.add(i)
        print('inequality ordering: {}'.format(self.inequality_ordering))

    def firstBasisWithBox(self):
        while not self.inside_box():
            for i, b in enumerate(self.B):
                if b not in self.boxIndices:
                    continue
                elif self.matrix[self.Row[i]][0] < 0:
                    # Primal infeasible Variable
                    self.i = i
                    break
            while self.matrix[self.Row[self.i]][self.Column[self.j]] <= 0:
                # To get primal feasible variable the sign of A[Row[i]][0] has to change
                # This happens if A[Row[i]][Column[j]] > 0
                self.j += 1
            if self.j == self.d:
                raise ValueError
            self.pivot()
        self.printInfo('After first basis with bounding box')

    def addBoxConstraints(self, constraints):
        """
        Dicts have to be initialized before.
        Constraints are of the form a0 + a1x1 + ... + a_d-1 x_d-1 >= 0
        """
        self.box_constraints = constraints
        self.matrix += constraints
        self.startBox = self.m + self.d
        self.boxIndices = set(list(range(self.startBox, self.startBox + len(constraints))))
        self.B += list(range(self.startBox, self.startBox + len(constraints)))
        self.C[-1] = self.m + self.d + len(constraints)
        self.Row += list(range(self.m + 1, self.m + 1 + len(constraints)))
        self.m += len(constraints)
        self.boxed = True
        self.inequality_ordering = list(range(self.d, self.d + self.m))

    def pivot_stays_in_box(self, i, j):
        pivotRow = self.Row[i]
        pivotColumn = self.Column[j]
        pivotElement = self.matrix[pivotRow][pivotColumn]
        insideBox = True
        for k, b in enumerate(self.B):
            if k == i:
                if self.C[j] not in self.boxIndices: # If not box variable is pivoted in we do not care about sign
                    print('Skipped because non pivot Variable', self.C[j])
                    continue
                elif self.computeEntryAfterPivot(self.Row[k], 0, pivotRow, pivotColumn, pivotElement) < 0:
                    insideBox = False
                    break

            if b in self.boxIndices and self.computeEntryAfterPivot(
                    self.Row[k], 0, pivotRow, pivotColumn, pivotElement
            ) < 0:
                insideBox = False
                break
        return insideBox

    @abstractmethod
    def select_pivot(self):
        pass

    def firstBoxBasisIndex(self):
        for i, b in enumerate(self.B):
            if b >= self.startBox:
                startBox = i
                return startBox

    def inside_box(self):
        for k, boxVar in enumerate(self.B):
            if boxVar not in self.boxIndices:
                continue
            if self.matrix[self.Row[k]][0] < 0:
                return False
        return True

    def search(self):
        self.printInfo('Search start:')
        self.i = self.d
        nextbasis = True
        backtrack = False
        while nextbasis:
            self.j = 0
            self.i = self.d
            while self.j < self.d or self.B[self.m] != self.m:
                if self.j == self.d - 1 and self.B[self.m] == self.m:
                    print('All bases found!')
                    print('bases:', self.bases)
                    print('vertices', self.vertices)
                    print('position vectors:', self.position_vectors)
                    nextbasis = False
                    break
                if backtrack:
                    print('Pivoting back!')
                    self.i, self.j = self.select_pivot()
                    self.pivot()
                    # self.i = self.d
                    # self.j += 1
                    self.increment()
                    print('i: {}, j: {}'.format(self.i, self.j))
                    backtrack = False
                else:
                    while self.j < self.d - 1 and not self.reverse():
                        self.increment()
                    if self.j == self.d - 1:
                        backtrack = True
                    else:
                        if self.lex_min():
                            self.appendSolution()
                        print('start tree search from new root')
                        break

    def appendSolution(self):
        print('Append basis: {}'.format(self.B))
        print('Vertex: {}'.format(self.getVertex()))
        self.bases.append(deepcopy(self.B))
        self.vertices.append(self.getVertex())
        self.position_vectors.append(self.getPositionVector())

    def getVertex(self):
        vertex = tuple(self.matrix[self.Row[k]][0] / self.det for k in range(1, self.d))
        return vertex

    def getPositionVector(self):
        basisIndex = self.d
        position_vector = []
        if self.boxed:
            last_real_constraint = self.startBox -1
        else:
            last_real_constraint = self.d + self.m -1
        for i in range(self.d, last_real_constraint + 1):
            if basisIndex <= self.m and self.B[basisIndex] == i:
                position_vector.append(-1 if self.matrix[self.Row[basisIndex]][0] < 0 else 1)
                basisIndex += 1
            else:
                position_vector.append(0)
        return position_vector

    def update(self):
        B_out = deepcopy(self.B[self.i])
        C_Out = deepcopy(self.C[self.j])
        self.B[self.i], self.C[self.j] = self.C[self.j], self.B[self.i]
        self.B, self.Row = self.sortDictionary(self.B, self.Row)
        self.C, self.Column = self.sortDictionary(self.C, self.Column)
        self.i = self.B.index(C_Out)
        self.j = self.C.index(B_out)

    def reverse(self):
        print('In reverse: i: {}, j:{}'.format(self.i, self.j))
        possibleReversePivot = self.necessaryConditionForReverse()
        if not possibleReversePivot:
            print('Not a possible reverse pivot!')
            return False
        self.pivot()
        i_forward, j_forward = self.select_pivot()
        if i_forward == self.i and j_forward == self.j:
            print('valid reverse')
            return True
        else:
            print('Not valid reverse: pivoting back')
            self.pivot()
            return False

    @abstractmethod
    def necessaryConditionForReverse(self):
        pass

    def lex_min(self):
        return True

    def printInfo(self, infoString=None):
        if infoString is not None:
            print(infoString)
        print('Basis: {}'.format(self.B))
        print('Cobasis: {}'.format(self.C))
        print('rows: {}, columns: {}'.format(self.Row, self.Column))
        print('det: {}'.format(self.det))
        print('matrix: ')
        self.pretty_print_matrix()

    def computeEntryAfterPivot(self, i, j, pivotRow, pivotColumn, pivotElement):
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
        print('pivot: outIndex: {}; inIndex: {}'.format(self.i, self.j))
        print('outVariable: {}; inVariable: {}'.format(self.B[self.i], self.C[self.j]))

        row = self.Row[self.i]
        column = self.Column[self.j]
        pivotElement = deepcopy(self.matrix[row][column])
        self.det = self.det if pivotElement > 0 else -self.det

        self.matrix = [
            [self.computeEntryAfterPivot(i, j, row, column, pivotElement) for j in range(self.d)]
            for i in range(self.m +1)
        ]
        self.det = pivotElement if pivotElement > 0 else - pivotElement
        self.update()
        self.printInfo('After pivot:')

    def sortDictionary(self, dictionary, locations):
        newDic, newLoc = zip(*sorted(zip(dictionary, locations)))
        return list(newDic), list(newLoc)

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

    @staticmethod
    def maxIndexSmallerNumber(list, number):
        """ Assumes sorted input list"""
        for i, element in enumerate(list):
            if element >= number:
                return i - 1