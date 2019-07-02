from reader import reader
from bland import Bland
from crisscross import CrissCross
from lrs import Lrs, SearchStatus
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

    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.variables_into_basis()


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
    lrs.augment_matrix_with_objective()
    lrs.init_dicts()
    lrs.add_box_constraints([boxConstraint1, boxConstraint2, boxConstraint3, boxConstraint4])
    lrs.variables_into_basis()
    lrs.move_into_box()
    lrs.set_objective()


def test_pivots_by_comparing():
    # Find vertices using bland
    b = Bland(*reader('data/cs_polytopes_boxed.ine'))
    b.augment_matrix_with_objective()
    b.init_dicts()
    b.add_box_constraints(b.bounding_box)
    b.first_basis()
    search = b.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()

    # Find vertices using criss-cross
    c = CrissCross(*reader('data/cs_polytopes_boxed.ine'))
    c.augment_matrix_with_objective()
    c.init_dicts()
    c.add_box_constraints(b.bounding_box)
    c.first_basis()
    search = c.search()
    status = SearchStatus.NONE
    while status != SearchStatus.DONE:
        status = search.__next__()

    assert len(set(b.vertices) - set(c.vertices)) == len(set(c.vertices) - set(b.vertices)) == 0

    from brute_force import brute_force_vertices

    vertices, vertices_inside_box = brute_force_vertices(b.hyperplanes, b.bounding_box)

    for vertex in set(b.vertices) - vertices_inside_box:
        assert any([cob.box_variable for cob in vertex.cobasis])
