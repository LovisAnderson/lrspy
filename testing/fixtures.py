from gmpy2 import mpz
from reader import reader
import pytest
from pathlib import Path

from lrs import Lrs


@pytest.fixture
def arrangement():
    bare_lrs = Lrs.__new__(Lrs)
    bare_lrs.matrix = [[mpz(0), mpz(0), mpz(0)],
                       [mpz(3), mpz(-1), mpz(-1)],
                       [mpz(-1), mpz(1), mpz(0)],
                       [mpz(-1), mpz(0), mpz(1)],
                       [mpz(-4), mpz(2), mpz(1)]]
    bare_lrs.B = [0, 3, 4, 5, 6]
    bare_lrs.C = [1, 2, 7]
    bare_lrs.Row = list(range(5))
    bare_lrs.Column = [1, 2, 0]
    bare_lrs.m = 4
    bare_lrs.d = 3
    bare_lrs.det = 1
    return bare_lrs


@pytest.fixture
def simplex():
    # 2d Simplex shifted by (1, 1)
    bare_lrs = Lrs.__new__(Lrs)
    bare_lrs.matrix = [[mpz(0), mpz(1), mpz(1)],
                       [mpz(-1), mpz(1), mpz(0)],
                       [mpz(-1), mpz(0), mpz(1)],
                       [mpz(2), mpz(-1), mpz(-1)], ]
    bare_lrs.B = [0, 3, 4, 5]
    bare_lrs.C = [1, 2, 6]
    bare_lrs.Row = list(range(4))
    bare_lrs.Column = [1, 2, 0]
    bare_lrs.m = 3
    bare_lrs.d = 3
    bare_lrs.det = 1

    return bare_lrs


@pytest.fixture
def lrs_from_file():
    p = Path(__file__).parents[1].joinpath('data/arrangement.ine')
    matrix, m, d = reader(str(p))
    lrs = Lrs(matrix, m, d)
    return lrs