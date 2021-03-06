from crisscross import CrissCross
from testing.fixtures import *
from lrs import SearchStatus
from brute_force import brute_force_vertices


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
    search = lrs.forest_search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    cobases = [[3, 4], [4, 5], [3, 5], [3, 6], [5, 6]]
    vertex_cobases = [vertex.cobasis for vertex in lrs.vertices]
    assert all(list(cobasis in vertex_cobases for cobasis in cobases))


def test_negative_search(from_file):
    lrs = CrissCross(*from_file)
    lrs.augment_matrix_with_objective()
    lrs.matrix[2][0] += 10
    lrs.init_dicts()
    lrs.first_basis()
    search = lrs.forest_search()
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
    search = lrs.forest_search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()


def test_box_search_from_file(from_file_boxed):
    lrs = CrissCross(*from_file_boxed)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.add_box_constraints(lrs.bounding_box)
    lrs.first_basis()
    search = lrs.forest_search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()


def test_zero_vertex(zero_vertex):
    lrs = CrissCross(*zero_vertex)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.first_basis()
    search = lrs.forest_search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    def is_zero(vertex):
        return all(vi == 0 for vi in vertex.coordinates)
    assert any(is_zero(v) for v in lrs.vertices)


def test_large_instance(nine_overlap):
    lrs = CrissCross(*nine_overlap)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.first_basis()
    search = lrs.forest_search()
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
    search = lrs.forest_search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    nr_vertices_box_interior = len([v for v in lrs.vertices
                                   if not any(c.box_variable for c in v.cobasis)])
    assert nr_vertices_box_interior == 70
    brute_vertices, brute_vertices_box_interior = brute_force_vertices(
        lrs.hyperplanes, lrs.box_constraints
    )
    assert nr_vertices_box_interior == len(brute_vertices_box_interior)


def test_hueh_boxed():
    from reader import reader
    path = '/nfs/OPTI/bzfander/lrs_python/data/hueh_358_boxed.ine'
    parsed = reader(path)
    lrs = CrissCross(*parsed)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.add_box_constraints(lrs.bounding_box)
    lrs.first_basis()
    search = lrs.forest_search()
    from lrs import SearchStatus
    status = SearchStatus.NONE
    i = 0
    while status != SearchStatus.DONE and i < 1000:
        status = search.__next__()
        i += 1
    assert i == 1000


def test_degenerated_2_boxed(arrangement_degenerated_2_boxed):
    lrs = CrissCross(*arrangement_degenerated_2_boxed)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.add_box_constraints(lrs.bounding_box)
    lrs.first_basis()
    search = lrs.forest_search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    nr_vertices_box_interior = len([v for v in lrs.vertices
                                    if not any(c.box_variable for c in v.cobasis)])
    brute_vertices, brute_vertices_box_interior = brute_force_vertices(
        lrs.hyperplanes, lrs.box_constraints
    )
    assert nr_vertices_box_interior == len(brute_vertices_box_interior)


def test_degenerated_2(arrangement_degenerated_2):
    lrs = CrissCross(*arrangement_degenerated_2)
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.first_basis()
    search = lrs.forest_search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()
    brute_vertices, _ = brute_force_vertices(
        lrs.hyperplanes, []
    )
    assert len(lrs.vertices) == len(brute_vertices)
