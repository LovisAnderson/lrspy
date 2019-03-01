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
        self.locations = [] # location[i] = (B/C, k) Variable i in Basis/Cobasis at index k
        self.m = m # Number of input hyperplanes
        self.d = d # Embedding dimension
        self.det = mpz(1) # determinant of the matrix, quasi the shared denominator
        # Initializing indices for search procedure
        self.bases = [] # list of bases found

    def augmentWithObjective(self):
        objectiveRow = [mpz(1)]*(self.d)
        objectiveRow[0] = mpz(0)
        self.matrix.insert(0, objectiveRow)

    def firstBasis(self):
        self.B = [0] + list(range(self.d , self.d + self.m))
        self.Row = list(range(self.m + 1))
        self.C = list(range(1, self.d)) + [self.d + self.m]
        self.Column = list(range(1, self.d)) + [0]

        for i in range(self.m+self.d + 1):
            if i in list(range(1, self.d )) + [self.d + self.m]:
                self.locations.append(('C', i))
            if i in [0] + list(range(self.d, self.d + self.m)):
                self.locations.append(('B', i))

        for inIndex in range(self.d - 1):
            outIndex = 1
            while (self.B[outIndex] in range(1, self.d) or
                   self.matrix[self.Row[outIndex]][self.Column[inIndex]] == 0):
                outIndex += 1
            self.pivot(outIndex, inIndex)

    def select_pivot(self):
        for i, dictIndicator in self.locations[1 : self.m + self.d]:
            if dictIndicator == 'B' and self.matrix[i][self.d] < 0:
                pass

    def search(self):
        i = 2
        j = 1
        while not (i > self.m and self.B[self.m] == self.m):
            while i <= self.m and not self.reverse():
                i, j = self.increment(i, j)
                if (i <= self.m):
                    if self.lex_min():
                        self.bases.append(self.B)
                    i = 2
                    j = 1
                else:
                    self.select_pivot()
                    self.pivot()
                    self.increment()

    def reverse(self, i, j):
        self.pivot(i, j)
        if self.select_pivot() == (j, i):
            return True
        else:
            self.pivot(i, j)
            return False

    def lex_min(self):
        pass

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
        self.locations[self.B[outIndex]] = ('B', outIndex)
        self.locations[self.C[inIndex]] = ('C', inIndex)
        self.printInfo('After pivot:')

    def sortDictionary(self, dictionary, locations):
        newDic, newLoc = zip(*sorted(zip(dictionary, locations)))
        return newDic, newLoc

    def increment(self, i, j):
        j += 1
        if j == self.n - self.m:
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
    from reader import reader, addSlacks

    matrix, m, d = reader('data/arrangement.ine')
    lrs = Lrs(matrix, m, d)
    lrs.augmentWithObjective()
    lrs.firstBasis()
    assert False

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
    bare_lrs.locations = [('B', 0), ('C', 0), ('C', 1), ('B', 1), ('B', 2), ('B', 3), ('B', 4),
                          ('C', 2)]
    return bare_lrs

def test_pivot(arrangement):
    from copy import deepcopy
    matrix_before = deepcopy(arrangement.matrix)
    B_before = deepcopy(arrangement.B)
    C_before = deepcopy(arrangement.C)
    arrangement.pivot(4, 0)
    assert arrangement.B == [0, 3, 4, 5, 1]
    assert arrangement.C == [6, 2, 7]
    arrangement.pivot(4,0)
    assert  matrix_before == arrangement.matrix
    assert B_before == arrangement.B
    assert C_before == arrangement.C
