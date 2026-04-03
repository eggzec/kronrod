# kronrod

**Gauss-Kronrod Quadrature Rules for Python**

[![Tests](https://github.com/eggzec/kronrod/actions/workflows/test.yml/badge.svg)](https://github.com/eggzec/kronrod/actions/workflows/test.yml)
[![Documentation](https://github.com/eggzec/kronrod/actions/workflows/docs.yml/badge.svg)](https://github.com/eggzec/kronrod/actions/workflows/docs.yml)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

[![codecov](https://codecov.io/github/eggzec/kronrod/graph/badge.svg)](https://codecov.io/github/eggzec/kronrod)
[![License: LGPL-2.1](https://img.shields.io/badge/License-LGPL--2.1-blue.svg)](LICENSE)

[![PyPI Downloads](https://img.shields.io/pypi/dm/kronrod.svg?label=PyPI%20downloads)](https://pypi.org/project/kronrod/)
[![Python versions](https://img.shields.io/pypi/pyversions/kronrod.svg)](https://pypi.org/project/kronrod/)

`kronrod` is a Python library for computing the abscissas and weights of
Gauss-Kronrod quadrature rules.  Given an N-point Gauss quadrature rule,
the library optimally adds N+1 points to produce a 2N+1 point
Gauss-Kronrod rule for numerical integration.

The advantage of using a Gauss and Gauss-Kronrod pair is that the second
rule, which uses 2*N+1 points, actually includes the N points in the
previous Gauss rule.  This means that the function values from that
computation can be reused.  This efficiency comes at the cost of a mild
reduction in the degree of polynomial precision of the Gauss-Kronrod rule.

## Quick example

```python
import kronrod

# Compute a 3-point Gauss / 7-point Gauss-Kronrod rule
x, w1, w2 = kronrod.kronrod(3)

print("Abscissas:", x)
print("Gauss-Kronrod weights:", w1)
print("Gauss weights:", w2)

# Adjust the rule from [-1, +1] to [0, 1]
x_adj, w1_adj, w2_adj = kronrod.kronrod_adjust(0.0, 1.0, 3, x, w1, w2)
```

## Installation

```bash
pip install kronrod
```

Requires Python 3.10+ and NumPy. See the
[full installation guide](https://eggzec.github.io/kronrod/installation/) for
uv, poetry, and source builds.

## Documentation

- [Theory](https://eggzec.github.io/kronrod/theory/) — Gauss-Kronrod quadrature, Kronrod extension, error estimation
- [Quickstart](https://eggzec.github.io/kronrod/quickstart/) — runnable examples
- [API Reference](https://eggzec.github.io/kronrod/api/) — function signatures and parameters
- [References](https://eggzec.github.io/kronrod/references/) — literature citations

## Routines

| Function | Description |
|---|---|
| `abwe1` | Calculate a Kronrod abscissa and weight |
| `abwe2` | Calculate a Gaussian abscissa and two weights |
| `kronrod` | Add N+1 points to an N-point Gaussian rule |
| `r8_abs` | Return the absolute value of a double |
| `timestamp` | Print the current YMDHMS date as a time stamp |

## References

- Robert Piessens, Maria Branders, "A Note on the Optimal Addition of Abscissas to Quadrature Formulas of Gauss and Lobatto," *Mathematics of Computation*, Volume 28, Number 125, January 1974, pages 135-139.

## Attribution

The original FORTRAN77 code was written by Robert Piessens and Maria Branders.
The C translation was made by
[John Burkardt](https://people.math.sc.edu/Burkardt/c_src/kronrod/kronrod.html)
and is distributed under the LGPL-2.1 license.

## License

LGPL-2.1 — see [LICENSE](LICENSE).
