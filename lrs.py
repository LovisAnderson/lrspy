from gmpy2 import mpz, divexact
import pytest
from copy import deepcopy


class Lrs:
    def __init__(self, inequality_matrix, m, d):
        self.matrix = inequality_matrix
        self.B = [] # Basis
        self.C = [] # Cobasis
        self.Row = [] # B[i] is to be found in Row[i] of matrix
        self.Column = [] # C[i] is to be found in Column[i] of matrix
        self.m = m # Number of input hyperplanes
        self.d = d # dimension of embedding space + 1
        self.det = mpz(1) # determinant of the matrix, quasi the shared denominator
        self.bases = [] # list of bases found
        self.vertices = []
        self.i = self.d
        self.j = 0

    def setObjective(self):
        for i in range(1, self.d):
            self.matrix[0][i] = mpz(-self.det)

    def augmentWithObjective(self):
        objectiveRow = [mpz(1)]*(self.d)
        objectiveRow[0] = mpz(0)
        self.matrix.insert(0, objectiveRow)

    def firstBasis(self):
        self.B = [0] + list(range(self.d , self.d + self.m))
        self.Row = list(range(self.m + 1))
        self.C = list(range(1, self.d)) + [self.d + self.m]
        self.Column = list(range(1, self.d)) + [0]

        for k in range(self.d - 1):
            self.j = 0
            self.i = 1
            while (self.B[self.i] in range(1, self.d) or
                   self.matrix[self.Row[self.i]][self.Column[self.j]] == 0):
                self.i += 1
            self.pivot()
        self.printInfo('After first basis')
        self.appendSolution()

    def select_pivot(self):
        basis_index = self.d
        cobasis_index = 0
        for i in range(self.d, self.d + self.m):
            if basis_index <= self.m and self.B[basis_index] == i:
                if self.matrix[self.Row[basis_index]][0] < 0:
                    # Bi is primal infeasible
                    print('B[{}] = {} primal infeasible'.format(basis_index, self.B[basis_index]))
                    for cobasis_index, c in enumerate(self.C):
                        if self.matrix[self.Row[basis_index]][self.Column[cobasis_index]] > 0:
                            return basis_index, cobasis_index
                    raise ValueError
                basis_index += 1
            elif cobasis_index < self.d and self.C[cobasis_index] == i:
                if self.matrix[0][self.Column[cobasis_index]] > 0:
                    # Ci is dual infeasible
                    print('C[{}] = {} dual infeasible'.format(cobasis_index, self.C[cobasis_index]))
                    for basis_index, b in enumerate(self.B):
                        if basis_index < self.d:
                            continue
                        if self.matrix[self.Row[basis_index]][self.Column[cobasis_index]] < 0:
                            return basis_index, cobasis_index
                    raise ValueError
                cobasis_index += 1

    def search(self):
        self.printInfo('Search start:')
        self.i = self.d
        nextbasis = True
        backtrack = False
        while nextbasis:
            self.j = 0
            while self.j < self.d or self.B[self.m] != self.m:
                if self.j == self.d - 1 and self.B[self.m] == self.m:
                    print('All bases found!')
                    print('bases:', self.bases)
                    print('vertices', self.vertices)
                    nextbasis = False
                    break
                if backtrack:
                    print('Pivoting back!')
                    self.i, self.j = self.select_pivot()
                    self.pivot()
                    self.increment()
                    backtrack = False
                else:
                    while self.j < self.d - 1 and not self.reverse:
                        self.increment()
                        print('Incrementing -> i={}, j={}'.format(self.i, self.j))
                    if self.j == self.d - 1:
                        backtrack = True
                    else:
                        if self.lex_min():
                            self.appendSolution()
                        print('start tree search from new root')
                        break

    def appendSolution(self):
        print('Append basis: {}'.format(self.B))
        self.bases.append(deepcopy(self.B))
        self.vertices.append(self.getVertex())

    def getVertex(self):
        vertex = tuple(self.matrix[self.Row[k]][0] / self.det for k in range(1, self.d))
        return vertex

    def update(self):
        B_out = deepcopy(self.B[self.i])
        C_Out = deepcopy(self.C[self.j])
        self.B[self.i], self.C[self.j] = self.C[self.j], self.B[self.i]
        self.B, self.Row = self.sortDictionary(self.B, self.Row)
        self.C, self.Column = self.sortDictionary(self.C, self.Column)
        self.i = self.B.index(C_Out)
        self.j = self.C.index(B_out)


    @property
    def reverse(self):
        print('In reverse: i: {}, j:{}'.format(self.i, self.j))
        possibleReversePivot = self.necessaryConditionForReverse()
        if not possibleReversePivot:
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


    def necessaryConditionForReverse(self):
        if self.matrix[self.Row[self.i]][0] > 0:
            if (
                    self.matrix[self.Row[self.i]][self.Column[self.j]] > 0 and
                    all(
                        self.matrix[self.Row[self.i]][self.Column[k]] <= 0
                        for k in range(0, maxIndexSmallerNumber(self.C, self.B[self.i]) + 1)
                    )
            ):
                return  True
        if self.matrix[0][self.Column[self.j]] < 0:
            if (self.matrix[self.Row[self.i]][self.Column[self.j]] < 0 and
                all(
                    self.matrix[self.Row[k]][self.Column[self.j]] <= 0
                    for k in range(1, maxIndexSmallerNumber(self.C, self.C[self.j]) + 1)
                )
            ):
                return True
        return False


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

def maxIndexSmallerNumber(list, number):
    """ Assumes sorted input list"""
    for i, element in enumerate(list):
        if element >= number:
            return i-1

def test_augment_with_objective():
    from reader import reader

    matrix, m, d = reader('data/arrangement.ine')
    lrs = Lrs(matrix, m, d)
    lrs.augmentWithObjective()
    lrs.firstBasis()

@pytest.fixture
def arrangement():
    bare_lrs = Lrs.__new__(Lrs)
    bare_lrs.matrix = [[mpz(0), mpz(0), mpz(0)],
                       [mpz(3), mpz(-1), mpz(-1)],
                       [mpz(-1), mpz(1), mpz(0)],
                       [mpz(-1), mpz(0), mpz(1)],
                       [mpz(-4), mpz(2), mpz(1)]]
    bare_lrs.B = [0, 3, 4, 5, 6]
    bare_lrs.C = [1, 2, 7]
    bare_lrs.Row = list(range(5))
    bare_lrs.Column = [1, 2, 0]
    bare_lrs.m = 4
    bare_lrs.d = 3
    bare_lrs.det = 1
    return bare_lrs

@pytest.fixture
def simplex():
    # 2d Simplex shifted by (1, 1)
    bare_lrs = Lrs.__new__(Lrs)
    bare_lrs.matrix = [[mpz(0), mpz(1), mpz(1)],
                       [mpz(-1), mpz(1), mpz(0)],
                       [mpz(-1), mpz(0), mpz(1)],
                       [mpz(2), mpz(-1), mpz(-1)],]
    bare_lrs.B = [0, 3, 4, 5]
    bare_lrs.C = [1, 2, 6]
    bare_lrs.Row = list(range(4))
    bare_lrs.Column = [1, 2, 0]
    bare_lrs.m = 3
    bare_lrs.d = 3
    bare_lrs.det = 1

    return bare_lrs


def test_pivot(arrangement):
    from copy import deepcopy
    matrix_before = deepcopy(arrangement.matrix)
    B_before = deepcopy(arrangement.B)
    C_before = deepcopy(arrangement.C)
    arrangement.i = 4
    arrangement.j = 0
    arrangement.pivot()
    assert list(arrangement.B) == [0, 1, 3, 4, 5]
    assert list(arrangement.C) == [2, 6, 7]
    # We pivot back and test if we get the same result
    arrangement.pivot()
    assert  matrix_before == arrangement.matrix
    assert B_before == arrangement.B
    assert C_before == arrangement.C


def test_select_pivot(simplex):
    bare_lrs = Lrs.__new__(Lrs)
    bare_lrs.matrix = [[mpz(0), mpz(1), mpz(1)],
                       [mpz(-1), mpz(1), mpz(0)],
                       [mpz(-1), mpz(0), mpz(1)],
                       [mpz(2), mpz(-1), mpz(-1)]]
    bare_lrs.B = [0, 1, 2, 5]
    bare_lrs.C = [3, 4, 6]
    bare_lrs.Row = list(range(4))
    bare_lrs.Column = [1, 2, 0]
    bare_lrs.m = 3
    bare_lrs.d = 3
    i, j = bare_lrs.select_pivot()
    assert i == 3
    assert j == 0

def test_search():
    from reader import reader

    matrix, m, d = reader('data/arrangement.ine')
    lrs = Lrs(matrix, m, d)
    lrs.augmentWithObjective()
    lrs.firstBasis()
    lrs.setObjective()
    lrs.search()
    assert False