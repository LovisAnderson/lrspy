from gmpy2 import mpz
from reader import reader
import pytest
from pathlib import Path
from lrs_datastructures import LrsDict, Variable

@pytest.fixture
def arrangement():
    matrix = [[mpz(0), mpz(0), mpz(0)],
              [mpz(3), mpz(-1), mpz(-1)],
              [mpz(-1), mpz(1), mpz(0)],
              [mpz(-1), mpz(0), mpz(1)],
              [mpz(-4), mpz(2), mpz(1)]]
    B = [0, 3, 4, 5, 6]
    B = LrsDict(list(Variable(b) for b in B))
    C = [1, 2, 7]
    C = LrsDict(list(Variable(c) for c in C))
    C.order = [1, 2, 0]
    m = 4
    d = 3
    det = 1
    attributes = {
        'matrix': matrix,
        'B': B,
        'C': C,
        'm': m,
        'd': d,
        'det': det,
        'nr_hyperplanes': m,
    }
    return attributes


@pytest.fixture
def arrangement2():
    matrix = [[mpz(2), mpz(1), mpz(0)],
              [mpz(2), mpz(-1), mpz(-1)],
              [mpz(1), mpz(1), mpz(0)],
              [mpz(-1), mpz(2), mpz(2)],
              [mpz(1), mpz(-1), mpz(-1)]]
    B = [0, 1, 2, 4, 6]
    B = LrsDict(list(Variable(b) for b in B))
    B[3].slack_variable = True
    B[3].hyperplane_index = 1
    B[4].slack_variable = True
    B[4].hyperplane_index = 3
    C = [3, 5, 7]
    C = LrsDict(list(Variable(c) for c in C))
    C.order = [2, 1, 0]
    C[0].slack_variable = True
    C[0].hyperplane_index = 0
    C[1].slack_variable = True
    C[1].hyperplane_index = 2
    d = 3
    m = 4
    boxed = False
    attributes = {
        'matrix': matrix,
        'B': B,
        'C': C,
        'm': m,
        'd': d,
        'boxed': boxed,
        'nr_hyperplanes': m,
    }

    return attributes

@pytest.fixture
def simplex():

    matrix = [[mpz(0), mpz(1), mpz(1)],
              [mpz(-1), mpz(1), mpz(0)],
              [mpz(-1), mpz(0), mpz(1)],
              [mpz(2), mpz(-1), mpz(-1)], ]
    B = [0, 3, 4, 5]
    B = LrsDict(list(Variable(b) for b in B))
    C = [1, 2, 6]
    C = LrsDict(list(Variable(c) for c in C))
    C.order = [1, 2, 0]
    m = 3
    d = 3
    det = 1

    attributes ={
        'matrix': matrix,
        'B': B,
        'C': C,
        'm': m,
        'd': d,
        'det': det,
        'nr_hyperplanes': m,
    }

    return attributes


@pytest.fixture
def from_file():
    p = Path(__file__).parents[1].joinpath('data/arrangement.ine')
    matrix, m, d = reader(str(p))
    return matrix, m, d