"""Kronrod — Gauss-Kronrod quadrature rule computation.

Compute the abscissas and weights of the 2N+1 point Gauss-Kronrod
quadrature formula obtained by optimally adding N+1 points to an
N-point Gauss quadrature rule.  The advantage of using a Gauss and
Gauss-Kronrod pair is that the second rule, which uses 2*N+1 points,
actually includes the N points in the previous Gauss rule so that the
function values from that computation can be reused.

Functions
---------
* :func:`kronrod` — add N+1 points to an N-point Gaussian rule
* :func:`kronrod_adjust` — adjust a rule from [-1, +1] to [a, b]

Example
-------
>>> import kronrod
>>> x, w1, w2 = kronrod.kronrod(3)
"""

from typing import TYPE_CHECKING

import numpy as np

from ._kronrod import kronrod as _kronrod


if TYPE_CHECKING:
    from numpy.typing import NDArray

_DEFAULT_EPS = 1.0e-6


def kronrod(
    n: int, eps: float = _DEFAULT_EPS
) -> tuple["NDArray[np.float64]", "NDArray[np.float64]", "NDArray[np.float64]"]:
    """Add N+1 Kronrod points to an N-point Gaussian rule.

    Compute the abscissas and weights of the 2N+1 point Gauss-Kronrod
    quadrature formula for integration over [-1, +1].  Only the
    non-negative abscissas are returned (the rule is symmetric).

    Parameters
    ----------
    n : int
        Order of the Gauss rule.
    eps : float, optional
        Requested absolute accuracy of the abscissas (default 1e-6).

    Returns
    -------
    x : numpy.ndarray
        Non-negative abscissas in decreasing order, shape ``(n+1,)``.
    w1 : numpy.ndarray
        Gauss-Kronrod weights, shape ``(n+1,)``.
    w2 : numpy.ndarray
        Gauss weights, shape ``(n+1,)`` (zero for Kronrod-only points).
    """
    x, w1, w2 = _kronrod(n, eps)
    return (
        np.asarray(x, dtype=np.float64),
        np.asarray(w1, dtype=np.float64),
        np.asarray(w2, dtype=np.float64),
    )


def kronrod_adjust(  # noqa: PLR0913, PLR0917
    a: float,
    b: float,
    n: int,
    x: "NDArray[np.float64]",
    w1: "NDArray[np.float64]",
    w2: "NDArray[np.float64]",
) -> tuple["NDArray[np.float64]", "NDArray[np.float64]", "NDArray[np.float64]"]:
    """Adjust a Gauss-Kronrod rule from [-1, +1] to [a, b].

    Parameters
    ----------
    a : float
        Left endpoint of the new interval.
    b : float
        Right endpoint of the new interval.
    n : int
        Order of the Gauss rule.
    x : numpy.ndarray
        Abscissas, shape ``(n+1,)``.
    w1 : numpy.ndarray
        Gauss-Kronrod weights, shape ``(n+1,)``.
    w2 : numpy.ndarray
        Gauss weights, shape ``(n+1,)``.

    Returns
    -------
    x : numpy.ndarray
        Adjusted abscissas, shape ``(n+1,)``.
    w1 : numpy.ndarray
        Adjusted Gauss-Kronrod weights, shape ``(n+1,)``.
    w2 : numpy.ndarray
        Adjusted Gauss weights, shape ``(n+1,)``.
    """
    x = np.array(x, dtype=np.float64, copy=True)
    w1 = np.array(w1, dtype=np.float64, copy=True)
    w2 = np.array(w2, dtype=np.float64, copy=True)
    half = (b - a) / 2.0
    x[: n + 1] = ((1.0 - x[: n + 1]) * a + (1.0 + x[: n + 1]) * b) / 2.0
    w1[: n + 1] = half * w1[: n + 1]
    w2[: n + 1] = half * w2[: n + 1]
    return x, w1, w2


__all__ = ["kronrod", "kronrod_adjust"]
