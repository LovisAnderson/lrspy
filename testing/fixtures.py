from gmpy2 import mpz
from reader import reader
import pytest
from pathlib import Path


@pytest.fixture
def arrangement():
    matrix = [[mpz(0), mpz(0), mpz(0)],
              [mpz(3), mpz(-1), mpz(-1)],
              [mpz(-1), mpz(1), mpz(0)],
              [mpz(-1), mpz(0), mpz(1)],
              [mpz(-4), mpz(2), mpz(1)]]
    B = [0, 3, 4, 5, 6]
    C = [1, 2, 7]
    Row = list(range(5))
    Column = [1, 2, 0]
    m = 4
    d = 3
    det = 1
    attributes = {
        'matrix': matrix,
        'B': B,
        'C': C,
        'Row': Row,
        'Column': Column,
        'm': m,
        'd': d,
        'det': det
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
    C = [3, 5, 7]
    Row = [0, 1, 2, 4, 3]
    Column = [2, 1, 0]
    d = 3
    m = 4
    boxed = False
    attributes = {
        'matrix': matrix,
        'B': B,
        'C': C,
        'Row': Row,
        'Column': Column,
        'm': m,
        'd': d,
        'boxed': boxed
    }

    return attributes

@pytest.fixture
def simplex():

    matrix = [[mpz(0), mpz(1), mpz(1)],
              [mpz(-1), mpz(1), mpz(0)],
              [mpz(-1), mpz(0), mpz(1)],
              [mpz(2), mpz(-1), mpz(-1)], ]
    B = [0, 3, 4, 5]
    C = [1, 2, 6]
    Row = list(range(4))
    Column = [1, 2, 0]
    m = 3
    d = 3
    det = 1

    attributes ={
        'matrix': matrix,
        'B': B,
        'C': C,
        'Row': Row,
        'Column': Column,
        'm': m,
        'd': d,
        'det': det
    }

    return attributes


@pytest.fixture
def from_file():
    p = Path(__file__).parents[1].joinpath('data/arrangement.ine')
    matrix, m, d = reader(str(p))
    return matrix, m, d