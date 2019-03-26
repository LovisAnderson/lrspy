from testing.fixtures import *
from gmpy2 import mpz
from bland import Bland
from lrs_datastructures import LrsDict, Variable


def test_select_pivot(simplex):
    bare_lrs = Bland.__new__(Bland)
    bare_lrs.matrix = [[mpz(0), mpz(1), mpz(1)],
                       [mpz(-1), mpz(1), mpz(0)],
                       [mpz(-1), mpz(0), mpz(1)],
                       [mpz(2), mpz(-1), mpz(-1)],
                       [mpz(1), mpz(-1), mpz(-1)]]
    bare_lrs.B = LrsDict([0, 1, 2, 5, 6])
    bare_lrs.C = LrsDict([3, 4, 7])
    bare_lrs.C.order = [1, 2, 0]
    bare_lrs.m = 4
    bare_lrs.d = 3
    bare_lrs.boxed = False
    i, j = bare_lrs.select_pivot()
    assert i == 4
    assert j == 0


def test_search(from_file):
    lrs = Bland(*from_file)
    lrs.augmentWithObjective()
    lrs.initDicts()
    lrs.firstBasis()
    lrs.setObjective()
    lrs.search()
    bases = [[0, 1, 2, 5, 6], [0, 1, 2, 3, 6], [0, 1, 2, 4, 6], [0, 1, 2, 4, 5], [0, 1, 2, 3, 4]]
    assert all(basis in lrs.bases for basis in bases)