from crisscross import CrissCross
from lrs_datastructures import LrsDict, Variable
from gmpy2 import mpz, mpfr
from testing.fixtures import *
from lrs import SearchStatus


def test_select_pivot(simplex):
    bare_lrs = CrissCross.__new__(CrissCross)
    bare_lrs.matrix = [[mpz(0), mpz(1), mpz(1)],
                       [mpz(-1), mpz(1), mpz(0)],
                       [mpz(-1), mpz(0), mpz(1)],
                       [mpz(2), mpz(-1), mpz(-1)]]

    bare_lrs.B = LrsDict([0, 1, 2, 5])
    bare_lrs.C = LrsDict([3, 4, 6])
    bare_lrs.C.order = [1, 2, 0]
    bare_lrs.m = 3
    bare_lrs.d = 3
    bare_lrs.boxed = False
    i, j = bare_lrs.select_pivot()
    assert i == 3
    assert j == 0


def test_search(from_file):
    lrs = CrissCross(*from_file)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.first_basis()
    lrs.set_objective()
    search = lrs.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    bases = [[0, 1, 2, 5, 6], [0, 1, 2, 3, 6], [0, 1, 2, 4, 6], [0, 1, 2, 4, 5], [0, 1, 2, 3, 4]]
    assert all(basis in lrs.bases for basis in bases)


def test_negative_search(from_file):
    lrs = CrissCross(*from_file)
    lrs.augment_matrix_with_objective()
    lrs.matrix[2][0] += 10
    lrs.info_string()
    lrs.init_dicts()
    lrs.first_basis()
    lrs.set_objective()
    lrs.search()


def test_box_search(from_file):
    lrs = CrissCross(*from_file)
    boxConstraint1 = [mpz(-5), mpz(4), mpz(0)]
    boxConstraint2 = [mpz(3), mpz(0), mpz(-1)]
    boxConstraint3 = [mpz(-1), mpz(0), mpz(2)]
    boxConstraint4 = [mpz(3), mpz(-1), mpz(0)]
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.add_box_constraints([boxConstraint1, boxConstraint2, boxConstraint3, boxConstraint4])
    lrs.first_basis()
    lrs.first_basis_with_box()
    lrs.set_objective()
    lrs.search()

def test_zero_vertex(zero_vertex):
    lrs = CrissCross(*zero_vertex)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.first_basis()
    lrs.set_objective()
    lrs.search()
    def is_zero(point):
        return all(vi == 0 for vi in point)
    assert any(is_zero(v) for v in lrs.vertices)
