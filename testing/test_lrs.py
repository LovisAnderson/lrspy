from reader import reader
from gmpy2 import mpz
from lrs import Lrs
from copy import deepcopy
from testing.fixtures import *


# Concrete dummy of abstract base class for testing purposes
class ConcreteLrs(Lrs):
    def __init__(self, hyperplane_matrix, m, d):
        super().__init__(hyperplane_matrix, m, d)

    def set_attributes(self, attributes):
        # Helper Method to set for test needed attributes
        for key, attr in attributes.items():
            setattr(self, key, attr)
ConcreteLrs.__abstractmethods__ = frozenset()


def test_augment_with_objective(from_file):

    matrix, m, d = from_file
    lrs = ConcreteLrs(matrix, m, d)

    lrs.augmentWithObjective()
    lrs.initDicts()
    lrs.firstBasis()


def test_position_vector(arrangement2):
    simpleLrs = ConcreteLrs.__new__(ConcreteLrs)
    simpleLrs.set_attributes(arrangement2)
    simpleLrs.getPositionVector()
    assert simpleLrs.getPositionVector() == [0, 1, 0, -1]


def test_pivot(arrangement):
    lrs = ConcreteLrs.__new__(ConcreteLrs)
    lrs.set_attributes(arrangement)
    matrix_before = deepcopy(lrs.matrix)
    B_before = deepcopy(lrs.B)
    C_before = deepcopy(lrs.C)
    lrs.i = 4
    lrs.j = 0
    lrs.pivot()
    assert list(lrs.B) == [0, 1, 3, 4, 5]
    assert list(lrs.C) == [2, 6, 7]
    # We pivot back and test if we get the same result
    lrs.pivot()
    assert matrix_before == lrs.matrix
    assert B_before == lrs.B
    assert C_before == lrs.C


def test_bounding_box(from_file):
    lrs = ConcreteLrs(*from_file)
    boxConstraint1 = [mpz(-5), mpz(4), mpz(0)]
    boxConstraint2 = [mpz(3), mpz(0), mpz(-1)]
    boxConstraint3 = [mpz(-1), mpz(0), mpz(2)]
    boxConstraint4 = [mpz(3), mpz(-1), mpz(0)]
    lrs.augmentWithObjective()
    lrs.initDicts()
    lrs.addBoxConstraints([boxConstraint1, boxConstraint2, boxConstraint3, boxConstraint4])
    lrs.firstBasis()
    lrs.firstBasisWithBox()
    lrs.setObjective()
