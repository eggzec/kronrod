"""Pytest tests for the kronrod Python package.

Test data is derived from the tabulated values in the legacy C test
suite (bin/legacy/kronrod_test.c) and the reference output
(bin/legacy/kronrod_test.txt).
"""

import pytest

import kronrod


# ---------------------------------------------------------------------------
# Tolerance
# ---------------------------------------------------------------------------
ATOL = 1.0e-5


# ---- test N=3 (odd) -------------------------------------------------------

# Expected compact output from legacy test (N=3, 4 entries)
N3_X_EXPECTED = [0.960491, 0.774597, 0.434244, 0.000000]
N3_W1_EXPECTED = [0.104656, 0.268488, 0.401397, 0.450917]
N3_W2_EXPECTED = [0.000000, 0.555556, 0.000000, 0.888889]

# Exact Gauss-Kronrod weights for N=3
N3_WK_EXACT = [
    0.104656226026467265194,
    0.268488089868333440729,
    0.401397414775962222905,
    0.450916538658474142345,
    0.401397414775962222905,
    0.268488089868333440729,
    0.104656226026467265194,
]

# Exact Gauss weights for N=3
N3_WG_EXACT = [
    0.555555555555555555556,
    0.888888888888888888889,
    0.555555555555555555556,
]


def test_kronrod_n3_abscissas() -> None:
    """N=3 non-negative abscissas match legacy output."""
    x, _, _ = kronrod.kronrod(3)
    for i, expected in enumerate(N3_X_EXPECTED):
        assert x[i] == pytest.approx(expected, abs=ATOL)


def test_kronrod_n3_gk_weights() -> None:
    """N=3 Gauss-Kronrod weights match legacy output."""
    _, w1, _ = kronrod.kronrod(3)
    for i, expected in enumerate(N3_W1_EXPECTED):
        assert w1[i] == pytest.approx(expected, abs=ATOL)


def test_kronrod_n3_gauss_weights() -> None:
    """N=3 Gauss weights match legacy output."""
    _, _, w2 = kronrod.kronrod(3)
    for i, expected in enumerate(N3_W2_EXPECTED):
        assert w2[i] == pytest.approx(expected, abs=ATOL)


def test_kronrod_n3_exact_gk_weights() -> None:
    """N=3 Gauss-Kronrod weights match exact values."""
    _, w1, _ = kronrod.kronrod(3)
    n = 3
    for i in range(1, 2 * n + 2):
        if i <= n + 1:
            i2 = i
        else:
            i2 = 2 * (n + 1) - i
        assert w1[i2 - 1] == pytest.approx(N3_WK_EXACT[i - 1], abs=1.0e-10)


def test_kronrod_n3_exact_gauss_weights() -> None:
    """N=3 Gauss weights match exact values."""
    _, _, w2 = kronrod.kronrod(3)
    n = 3
    for i in range(1, n + 1):
        if 2 * i <= n + 1:
            i2 = 2 * i
        else:
            i2 = 2 * (n + 1) - 2 * i
        assert w2[i2 - 1] == pytest.approx(N3_WG_EXACT[i - 1], abs=1.0e-10)


# ---- test N=4 (even) ------------------------------------------------------

N4_X_EXPECTED = [0.976560, 0.861136, 0.640286, 0.339981, 0.000000]
N4_W1_EXPECTED = [0.062977, 0.170054, 0.266798, 0.326949, 0.346443]
N4_W2_EXPECTED = [0.000000, 0.347855, 0.000000, 0.652145, 0.000000]


def test_kronrod_n4_abscissas() -> None:
    """N=4 non-negative abscissas match legacy output."""
    x, _, _ = kronrod.kronrod(4)
    for i, expected in enumerate(N4_X_EXPECTED):
        assert x[i] == pytest.approx(expected, abs=ATOL)


def test_kronrod_n4_gk_weights() -> None:
    """N=4 Gauss-Kronrod weights match legacy output."""
    _, w1, _ = kronrod.kronrod(4)
    for i, expected in enumerate(N4_W1_EXPECTED):
        assert w1[i] == pytest.approx(expected, abs=ATOL)


def test_kronrod_n4_gauss_weights() -> None:
    """N=4 Gauss weights match legacy output."""
    _, _, w2 = kronrod.kronrod(4)
    for i, expected in enumerate(N4_W2_EXPECTED):
        assert w2[i] == pytest.approx(expected, abs=ATOL)


# ---- test integral estimate ------------------------------------------------


def _test_func(x: float) -> float:
    return 1.0 / (x * x + 1.005)


EXACT_INTEGRAL = 1.5643964440690497731
_MAX_ORDER = 25
_CONVERGENCE_TOL = 0.0001
_ACCURACY_TOL = 1.0e-4


def test_integral_convergence() -> None:
    """Gauss-Kronrod pair converges for 1/(x^2+1.005) on [-1,1]."""
    n = 1
    while n <= _MAX_ORDER:
        x, w1, w2 = kronrod.kronrod(n)

        i1 = float(w1[n] * _test_func(float(x[n])))
        i2 = float(w2[n] * _test_func(float(x[n])))

        for i in range(n):
            xi = float(x[i])
            i1 += float(w1[i]) * (_test_func(-xi) + _test_func(xi))
            i2 += float(w2[i]) * (_test_func(-xi) + _test_func(xi))

        if abs(i1 - i2) < _CONVERGENCE_TOL:
            assert abs(EXACT_INTEGRAL - i1) < _ACCURACY_TOL
            return

        n = 2 * n + 1

    pytest.fail(f"Integral did not converge within N={_MAX_ORDER}")


def test_integral_converges_at_n7() -> None:
    """Legacy output shows convergence at N=7."""
    n = 7
    x, w1, w2 = kronrod.kronrod(n)

    i1 = float(w1[n] * _test_func(float(x[n])))
    i2 = float(w2[n] * _test_func(float(x[n])))

    for i in range(n):
        xi = float(x[i])
        i1 += float(w1[i]) * (_test_func(-xi) + _test_func(xi))
        i2 += float(w2[i]) * (_test_func(-xi) + _test_func(xi))

    assert abs(i1 - i2) < _CONVERGENCE_TOL
    assert abs(EXACT_INTEGRAL - i1) < _ACCURACY_TOL


# ---- test kronrod_adjust --------------------------------------------------


def test_kronrod_adjust_interval() -> None:
    """Adjusted abscissas lie within [a, b]."""
    x, w1, w2 = kronrod.kronrod(3)
    x_adj, _, _ = kronrod.kronrod_adjust(0.0, 2.0, 3, x, w1, w2)
    adj_lo = -0.01
    adj_hi = 2.01
    for xi in x_adj:
        assert adj_lo <= float(xi) <= adj_hi


def test_kronrod_adjust_does_not_modify_original() -> None:
    """kronrod_adjust returns new arrays, does not modify originals."""
    x, w1, w2 = kronrod.kronrod(3)
    x_orig = x.copy()
    kronrod.kronrod_adjust(0.0, 1.0, 3, x, w1, w2)
    for i in range(len(x)):
        assert float(x[i]) == pytest.approx(float(x_orig[i]), abs=1e-15)


# ---- test output shapes ---------------------------------------------------


@pytest.mark.parametrize("n", [1, 3, 4, 7, 10])
def test_output_shapes(n: int) -> None:
    """Output arrays have shape (n+1,)."""
    x, w1, w2 = kronrod.kronrod(n)
    assert len(x) == n + 1
    assert len(w1) == n + 1
    assert len(w2) == n + 1
