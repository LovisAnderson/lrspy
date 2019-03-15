from crisscross import CrissCross
from gmpy2 import mpz
from testing.fixtures import *


def test_select_pivot(simplex):
    bare_lrs = CrissCross.__new__(CrissCross)
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
    bare_lrs.boxed = False
    i, j = bare_lrs.select_pivot()
    assert i == 3
    assert j == 0


def test_search(from_file):
    lrs = CrissCross(*from_file)
    lrs.augmentWithObjective()
    lrs.initDicts()
    lrs.firstBasis()
    lrs.setObjective()
    lrs.search()
    bases = [[0, 1, 2, 5, 6], [0, 1, 2, 3, 6], [0, 1, 2, 4, 6], [0, 1, 2, 4, 5], [0, 1, 2, 3, 4]]
    assert all(basis in lrs.bases for basis in bases)


def test_negative_search(from_file):
    lrs = CrissCross(*from_file)
    lrs.augmentWithObjective()
    lrs.matrix[2][0] += 10
    lrs.printInfo()
    lrs.initDicts()
    lrs.firstBasis()
    lrs.setObjective()
    lrs.search()