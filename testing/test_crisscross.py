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
    search = lrs.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    cobases = [[3, 4, 7], [4, 5, 7], [3, 5, 7], [3, 6, 7], [5, 6, 7]]
    assert all(cobasis in lrs.cobases for cobasis in cobases)


def test_negative_search(from_file):
    lrs = CrissCross(*from_file)
    lrs.augment_matrix_with_objective()
    lrs.matrix[2][0] += 10
    lrs.info_string()
    lrs.init_dicts()
    lrs.first_basis()
    search = lrs.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()


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
    search = lrs.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()


def test_box_search_from_file(from_file_boxed):
    lrs = CrissCross(*from_file_boxed)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.add_box_constraints(lrs.bounding_box)
    lrs.first_basis()
    search = lrs.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()


def test_zero_vertex(zero_vertex):
    lrs = CrissCross(*zero_vertex)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.first_basis()
    search = lrs.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    def is_zero(point):
        return all(vi == 0 for vi in point)
    assert any(is_zero(v) for v in lrs.vertices)


def test_large_instance(nine_overlap):
    lrs = CrissCross(*nine_overlap)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.first_basis()
    search = lrs.search()
    status = SearchStatus.NONE
    i = 0
    while status != SearchStatus.DONE and i <= 100:
        status = search.__next__()
        i += 1


def test_cs_boxed(cs_polytopes_boxed):
    lrs = CrissCross(*cs_polytopes_boxed)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.add_box_constraints(lrs.bounding_box)
    lrs.first_basis()
    search = lrs.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    nr_vertices_in_box = 0
    for cobasis in lrs.cobases:
        if not any(c.box_variable for c in cobasis[:-1]):
            nr_vertices_in_box += 1
    assert nr_vertices_in_box == 70


def test_hueh_boxed():
    from reader import reader
    path = '/nfs/OPTI/bzfander/lrs_python/data/hueh_358_boxed.ine'
    parsed = reader(path)
    lrs = CrissCross(*parsed)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.add_box_constraints(lrs.bounding_box)
    lrs.first_basis()
    search = lrs.search()
    from lrs import SearchStatus
    status = SearchStatus.NONE
    i = 0
    while status != SearchStatus.DONE and i < 1000:
        status = search.__next__()
        i += 1
    assert i == 1000