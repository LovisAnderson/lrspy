from testing.fixtures import *
from gmpy2 import mpz
from bland import Bland
from lrs_datastructures import LrsDict, Variable
from lrs import SearchStatus


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
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.first_basis()
    search = lrs.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()

    cobases = [[3, 4], [4, 5], [3, 5], [3, 6], [5, 6]]
    vertex_cobases = [vertex.cobasis for vertex in lrs.vertices]
    assert all(list(cobasis in vertex_cobases for cobasis in cobases))


def test_cs_boxed(cs_polytopes_boxed):
    lrs = Bland(*cs_polytopes_boxed)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.add_box_constraints(lrs.bounding_box)
    lrs.first_basis()
    search = lrs.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    nr_vertices_in_box = 0
    vertex_cobases = [vertex.cobasis for vertex in lrs.vertices]
    for cobasis in vertex_cobases:
        if not any(c.box_variable for c in cobasis):
            nr_vertices_in_box += 1
    assert nr_vertices_in_box == 70