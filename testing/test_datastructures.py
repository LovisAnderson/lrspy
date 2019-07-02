from lrs_datastructures import LrsDict, Variable
from gmpy2 import mpfr, mpz
from testing.fixtures import *
from vertex import Vertex

def test_dict_sort():
    a = LrsDict()
    a.append(1)
    a.append(4)
    a.append(2)
    a.order = [0, 1, 2]
    a.sort_respecting_order()
    assert a == [1, 2, 4]
    assert a.order == [0, 2, 1]


def test_position_vector(arrangement2):
    vertex = Vertex((mpfr(1), mpfr(2)))
    matrix = [
        [mpz(0), mpz(-1), mpz(-1)],
        [mpz(1), mpz(0), mpz(1)],
        [mpz(2), mpz(-1), mpz(-1)],
        [mpz(1), mpz(-1), mpz(-1)],
        [mpz(1), mpz(0), mpz(-2)]
    ]
    d = 3
    B = [0, 1, 2, 3, 4]
    B = [Variable(b) for b in B]
    for i, b in enumerate(B[d:]):
        b.slack_variable = True
        b.box_variable = False
        b.hyperplane_index = i
    B = LrsDict(B)
    B.order = [0, 1, 2, 3, 4]
    vertex.compute_position_vector(matrix, B, 2, d)
    assert vertex.position_vector == [1, 1]