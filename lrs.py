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

        for i in range(self.d - 1):
            inIndex = 0
            outIndex = 1
            while (self.B[outIndex] in range(1, self.d) or
                   self.matrix[self.Row[outIndex]][self.Column[inIndex]] == 0):
                outIndex += 1
            self.pivot(outIndex, inIndex)
        self.printInfo('After first basis')

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
                if self.matrix[0][self.Column[cobasis_index]] < 0:
                    # Ci is dual infeasible
                    print('C[{}] = {} dual infeasible'.format(cobasis_index, self.C[cobasis_index]))
                    for basis_index, b in enumerate(self.B):
                        if basis_index < self.d:
                            continue
                        if self.matrix[self.Row[basis_index][self.Column[cobasis_index]]] < 0:
                            return basis_index, cobasis_index
                    raise ValueError
                cobasis_index += 1

    def search(self):
        i = self.d
        j = 0
        while (j < self.d or self.B[self.m] != self.m):
            while i <= self.m and not self.reverse(i, j):
                i, j = self.increment(i, j)
            if (i <= self.m):
                if self.lex_min():
                    self.bases.append(self.B)
                i = self.d
                j = 0
            else:
                self.select_pivot()
                self.pivot()
                self.increment()

    def reverse(self, i, j):
        b = self.B[i]
        c = self.C[j]
        self.pivot(i, j)
        i_forward, j_forward = self.select_pivot()
        if self.B[i_forward] == c and self.C[j_forward] == b:
            return True
        else:
            i_back = self.B.index(c)
            j_back = self.C.index(b)
            self.pivot(i_back, j_back)
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

    def pivot(self, outIndex, inIndex):
        self.printInfo('Before pivot')
        print('pivot: outIndex: {}; inIndex: {}'.format(outIndex, inIndex))
        print('outVariable: {}; inVariable: {}'.format(self.B[outIndex], self.C[inIndex]))

        row = self.Row[outIndex]
        column = self.Column[inIndex]
        pivotElement = deepcopy(self.matrix[row][column])
        self.det = self.det if pivotElement > 0 else -self.det

        self.matrix = [
            [self.computeEntryAfterPivot(i, j, row, column, pivotElement) for j in range(self.d)]
            for i in range(self.m +1)
        ]
        self.det = pivotElement if pivotElement > 0 else - pivotElement
        self.B[outIndex], self.C[inIndex] = self.C[inIndex], self.B[outIndex]
        self.B, self.Row = self.sortDictionary(self.B, self.Row)
        self.C, self.Column = self.sortDictionary(self.C, self.Column)
        self.printInfo('After pivot:')

    def sortDictionary(self, dictionary, locations):
        newDic, newLoc = zip(*sorted(zip(dictionary, locations)))
        return list(newDic), list(newLoc)

    def increment(self, i, j):
        j += 1
        if j == self.d:
            j = 1
            i += 1
        return i, j

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


def test_pivot(arrangement):
    from copy import deepcopy
    matrix_before = deepcopy(arrangement.matrix)
    B_before = deepcopy(arrangement.B)
    C_before = deepcopy(arrangement.C)
    arrangement.pivot(4, 0)
    assert list(arrangement.B) == [0, 1, 3, 4, 5]
    assert list(arrangement.C) == [2, 6, 7]
    # We pivot back and test if we get the same result
    arrangement.pivot(1,1)
    assert  matrix_before == arrangement.matrix
    assert B_before == arrangement.B
    assert C_before == arrangement.C


def test_select_pivot(arrangement):
    i, j = arrangement.select_pivot()
    assert i == 2
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