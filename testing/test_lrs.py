from reader import reader
from gmpy2 import mpz
from lrs import Lrs
from copy import deepcopy
from testing.fixtures import *


def test_augment_with_objective():
    from reader import reader

    matrix, m, d = reader('data/arrangement.ine')
    lrs = Lrs(matrix, m, d)
    lrs.augmentWithObjective()
    lrs.initDicts()
    lrs.firstBasis()


def test_position_vector():
    simpleLrs = Lrs.__new__(Lrs)
    simpleLrs.matrix = [[mpz(2), mpz(1), mpz(0)],
                        [mpz(2), mpz(-1), mpz(-1)],
                        [mpz(1), mpz(1), mpz(0)],
                        [mpz(-1), mpz(2), mpz(2)],
                        [mpz(1), mpz(-1), mpz(-1)]]
    simpleLrs.B = [0, 1, 2, 4, 6]
    simpleLrs.C = [3, 5, 7]
    simpleLrs.Row = [0, 1, 2, 4, 3]
    simpleLrs.Column = [2, 1, 0]
    simpleLrs.d = 3
    simpleLrs.m = 4
    simpleLrs.getPositionVector()
    assert simpleLrs.getPositionVector() == [0, 1, 0, -1]


def test_pivot(arrangement):
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
    assert matrix_before == arrangement.matrix
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


def test_search(lrs_from_file):
    lrs = lrs_from_file
    lrs.augmentWithObjective()
    lrs.initDicts()
    lrs.firstBasis()
    lrs.setObjective()
    lrs.search()
    bases = [[0, 1, 2, 5, 6], [0, 1, 2, 3, 6], [0, 1, 2, 4, 6], [0, 1, 2, 4, 5], [0, 1, 2, 3, 4]]
    assert all(basis in lrs.bases for basis in bases)


def test_negative_search(lrs_from_file):
    lrs = lrs_from_file
    lrs.augmentWithObjective()
    lrs.matrix[2][0] += 10
    lrs.printInfo()
    lrs.initDicts()
    lrs.firstBasis()
    lrs.setObjective()
    lrs.search()

def test_bounding_box(lrs_from_file):
    lrs = lrs_from_file
    boxConstraint1 = [mpz(-5), mpz(4), mpz(0)]
    boxConstraint2 = [mpz(3), mpz(0), mpz(-1)]
    boxConstraint3 = [mpz(-1), mpz(0), mpz(2)]
    boxConstraint4 = [mpz(3), mpz(0), mpz(-1)]
    lrs.addBoxConstraints([boxConstraint1, boxConstraint2, boxConstraint3, boxConstraint4])
    print(lrs.printInfo())
    assert False