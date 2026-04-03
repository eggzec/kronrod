# API Reference

`kronrod` provides a Python interface to the Gauss-Kronrod quadrature C
library (Piessens & Branders, 1974; C translation by Burkardt) for computing
abscissas and weights of Gauss-Kronrod rules.  The numerical core is
implemented in C and compiled via `f2py`, providing near-native performance
with a clean Python interface.  See the [Theory](theory.md) and
[Quickstart](quickstart.md) for mathematical background and usage examples.

---

## Main Features

- Gauss-Kronrod rule computation: `kronrod`
- Interval adjustment: `kronrod_adjust`
- Embedded Gauss points reused in Kronrod extension
- Configurable accuracy parameter `eps`

## Function Summary

| Function | Category | Computes |
|---|---|---|
| `kronrod(n, eps)` | Rule computation | Abscissas and weights for 2N+1 point Gauss-Kronrod rule |
| `kronrod_adjust(a, b, n, x, w1, w2)` | Interval mapping | Adjust rule from $[-1, +1]$ to $[a, b]$ |

---

## `kronrod(n, eps=1e-6)`

Add N+1 Kronrod points to an N-point Gaussian rule.

Compute the abscissas and weights of the 2N+1 point Gauss-Kronrod
quadrature formula for integration over [-1, +1].  Only the non-negative
abscissas are returned (the rule is symmetric).

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | Order of the Gauss rule |
| `eps` | `float` | Requested absolute accuracy of the abscissas (default 1e-6) |
| **return** | `tuple[ndarray, ndarray, ndarray]` | `(x, w1, w2)` — abscissas, Gauss-Kronrod weights, Gauss weights |

```python
import kronrod

x, w1, w2 = kronrod.kronrod(3)
print(x)  # [0.96049127 0.77459667 0.43424375 0.        ]
print(w1)  # [0.10465623 0.26848809 0.40139741 0.45091654]
print(w2)  # [0.         0.55555556 0.         0.88888889]
```

---

## `kronrod_adjust(a, b, n, x, w1, w2)`

Adjust a Gauss-Kronrod rule from [-1, +1] to [a, b].

Linearly maps the abscissas and scales the weights so that the rule
integrates over the interval [a, b] instead of [-1, +1].

| Parameter | Type | Description |
|---|---|---|
| `a` | `float` | Left endpoint of the new interval |
| `b` | `float` | Right endpoint of the new interval |
| `n` | `int` | Order of the Gauss rule |
| `x` | `ndarray` | Abscissas from `kronrod()` |
| `w1` | `ndarray` | Gauss-Kronrod weights from `kronrod()` |
| `w2` | `ndarray` | Gauss weights from `kronrod()` |
| **return** | `tuple[ndarray, ndarray, ndarray]` | `(x, w1, w2)` — adjusted abscissas and weights |

```python
import kronrod

x, w1, w2 = kronrod.kronrod(3)
x_adj, w1_adj, w2_adj = kronrod.kronrod_adjust(0.0, 1.0, 3, x, w1, w2)
print(x_adj)  # adjusted abscissas in [0, 1]
```

---

## List of C routines

| Routine | Description |
|---|---|
| `abwe1` | Calculate a Kronrod abscissa and weight |
| `abwe2` | Calculate a Gaussian abscissa and two weights |
| `kronrod` | Add N+1 points to an N-point Gaussian rule |
| `r8_abs` | Return the absolute value of a double |
| `timestamp` | Print the current YMDHMS date as a time stamp |
