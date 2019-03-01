from enum import Enum
from gmpy2 import mpq, mpz


def reader(path):
    readMatrix = False
    readDimensions = False
    numberType = NumberType.undefined
    matrix = []
    with open(path) as infile:
        for line in infile.readlines():
            if line.startswith('begin'):
                readMatrix = True
                readDimensions = True
                continue
            elif line.startswith('end'):
                readMatrix = False
            elif readDimensions:
                m, d, numberType= parseMatrixMeta(line)
                readDimensions = False
                continue
            elif readMatrix:
                matrix.append(parseLine(line, d, numberType))
    return matrix, m, d


def addSlacks(matrix, m):
    for i, row in enumerate(matrix):
        slack = [mpq(0)]*m
        slack[i] = mpq(1)
        row += slack
    return matrix


def parseMatrixMeta(line):
    m, d, numberTypeStr = line.split()
    m = int(m)
    d = int(d)
    numberType = NumberType[numberTypeStr]
    return m, d, numberType


def parseLine(line, dim, numberType):
    entries = line.split()
    if len(entries) != dim:
        print(len(entries))
        raise ValueError('Line has to have {} many entries: {}!'.format(dim, entries))
    if numberType == NumberType.float:
        raise NotImplementedError
    if numberType == NumberType.rational:
        matrix_row = []
        for entry in entries:
            splitted = entry.split('/')
            if len(splitted) not in [1, 2]:
                raise ValueError
            if len(splitted) == 1:
                matrix_row.append(mpz(int(splitted[0])))
            else:
                matrix_row.append(mpq(int(splitted[0]), int(splitted[1])))
    if numberType == NumberType.integer:
        matrix_row = [int(a) for a in line]
    return matrix_row


class NumberType(Enum):
    undefined = 0
    rational = 1
    integer = 2
    float = 3

    def __str__(self):
        return self.name


def test_parse():
    A, m, d = reader('data/arrangement.ine')
    assert m == 4
    assert d == 3
    A1 = [
        [3, -1, -1],
        [-1, 1, 0],
        [-1, 0, 1],
        [3, -2, 0]
        ]
    A1 = [[mpq(a) for a in ai] for ai in A1]
    assert A == A1
    A = addSlacks(A, m)
    A1[0] += [mpz(1), mpz(0), mpz(0), mpz(0)]
    A1[1] += [mpz(0), mpz(1), mpz(0), mpz(0)]
    A1[2] += [mpz(0), mpz(0), mpz(1), mpz(0)]
    A1[3] += [mpz(0), mpz(0), mpz(0), mpz(1)]
    assert A == A1
