# kronrod

**Gauss-Kronrod Quadrature Rule Computation for Python**

---

## Overview

`kronrod` is a Python library for computing the abscissas and weights of
Gauss-Kronrod quadrature rules.  Given an N-point Gauss quadrature rule,
the library optimally adds N+1 points to produce a 2N+1 point
Gauss-Kronrod rule for numerical integration over [-1, +1].

The advantage of using a Gauss and Gauss-Kronrod pair is that the second
rule, which uses 2*N+1 points, actually includes the N points in the
previous Gauss rule.  This means that the function values from that
computation can be reused.  This efficiency comes at the cost of a mild
reduction in the degree of polynomial precision of the Gauss-Kronrod rule.

The algorithm was originally developed by Robert Piessens and Maria Branders
(1974) and is widely used in adaptive quadrature libraries such as QUADPACK.

## Requirements

- [NumPy](http://www.numpy.org/)

## Example Usage

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

## Main Features

1. **Gauss-Kronrod rule computation** — `kronrod` adds N+1 Kronrod abscissas to an N-point Gauss rule.
2. **Interval adjustment** — `kronrod_adjust` maps a rule from [-1, +1] to any [a, b].
3. **Embedded Gauss points** — the Kronrod extension reuses all N Gauss points, saving function evaluations.
4. Configurable accuracy parameter `eps` for abscissa precision.

## See Also

- [NumPy](http://www.numpy.org/): Array library providing the numerical foundation
- [SciPy](http://scipy.org): General scientific computing — includes adaptive quadrature routines

## Acknowledgements

The author thanks Robert Piessens, Maria Branders, and
[John Burkardt](https://people.math.sc.edu/Burkardt/c_src/kronrod/kronrod.html)
for the original FORTRAN77 and C implementations of the Gauss-Kronrod
algorithm, which provide the numerical core of this package.

## References

- Piessens, R. and Branders, M., "A Note on the Optimal Addition of Abscissas to Quadrature Formulas of Gauss and Lobatto," *Mathematics of Computation*, 28(125), 135–139, January 1974.
