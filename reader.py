from enum import Enum
from gmpy2 import mpq, mpz


def reader(path):
    readMatrix = False
    readDimensions = False
    offset_at_the_end = False
    numberType = NumberType.undefined
    matrix = []
    with open(path) as infile:
        for line in infile.readlines():
            if line.startswith('flipped representation'):
                offset_at_the_end = True
            if line.startswith('begin'):
                readMatrix = True
                readDimensions = True
                continue
            elif line.startswith('end'):
                readMatrix = False
            elif readDimensions:
                m, d, numberType= parse_matrix_meta(line)
                readDimensions = False
                continue
            elif readMatrix:
                matrix.append(parse_line(line, d, numberType, offset_at_the_end))
    return matrix, m, d


def addSlacks(matrix, m):
    for i, row in enumerate(matrix):
        slack = [mpq(0)]*m
        slack[i] = mpq(1)
        row += slack
    return matrix

def parse_matrix_meta(line):
    m, d, numberTypeStr = line.split()
    m = int(m)
    d = int(d)
    numberType = NumberType[numberTypeStr]
    return m, d, numberType


def parse_line(line, dim, numberType, offset_at_the_end=False):
    entries = line.split()
    if offset_at_the_end:
        entries = [entries[-1]] + entries[:-1]
    if len(entries) != dim:
        print(len(entries))
        raise ValueError('Line has to have {} many entries: {}!'.format(dim, entries))
    if numberType == NumberType.float:
        print(entries)
        matrix_row = [upscaled_mpz_from_float(float(entry)) for entry in entries]
        print(matrix_row)
    elif numberType == NumberType.rational:
        matrix_row = []
        for entry in entries:
            splitted = entry.split('/')
            if len(splitted) not in [1, 2]:
                raise ValueError
            if len(splitted) == 1:
                matrix_row.append(mpz(int(splitted[0])))
            else:
                matrix_row.append(mpq(int(splitted[0]), int(splitted[1])))
    elif numberType == NumberType.integer:
        matrix_row = [mpz(a) for a in line]

    return matrix_row

def upscaled_mpz_from_float(entry, scaling=1e+7):
    return mpz(entry*scaling)

class NumberType(Enum):
    undefined = 0
    rational = 1
    integer = 2
    float = 3

    def __str__(self):
        return self.name


def gcd(a, b):
    """Return greatest common divisor using Euclid's Algorithm."""
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    """Return lowest common multiple."""
    return a * b // gcd(a, b)


def lcmm(*args):
    """Return lcm of args."""
    z = lcm(args[0], args[1])
    if len(args) == 2:
        return z
    for arg in args[2:]:
        z = lcm(z, arg)
    return z

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
