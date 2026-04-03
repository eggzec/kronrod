# Theory

## Overview

`kronrod` implements the **Gauss-Kronrod** quadrature rule construction
algorithm of Piessens and Branders (1974).  Given an N-point Gauss
quadrature rule, the algorithm optimally adds N+1 points to produce a
2N+1 point Gauss-Kronrod rule for numerical integration.

---

## Background: Numerical integration

Numerical integration (quadrature) approximates a definite integral by a
weighted sum of function values:

$$
\int_{a}^{b} f(x)\,dx \approx \sum_{i=1}^{m} w_i\, f(x_i).
$$

Classical approaches include:

- **Newton-Cotes rules** — equally spaced nodes, low order
- **Gauss quadrature** — optimal node placement, exact for polynomials of degree $2N-1$
- **Gauss-Kronrod** — extends a Gauss rule by reusing its nodes, enabling efficient error estimation

Gauss-Kronrod rules are the foundation of adaptive quadrature libraries
such as QUADPACK.

---

## Gauss quadrature

An N-point Gauss quadrature rule on $[-1, +1]$ chooses nodes $x_i$ and
weights $w_i$ so that:

$$
\int_{-1}^{1} f(x)\,dx \approx \sum_{i=1}^{N} w_i\, f(x_i)
$$

is exact for all polynomials of degree $\le 2N - 1$.  The nodes are the
zeros of the Legendre polynomial $P_N(x)$.

---

## The Kronrod extension

### Motivation

In adaptive quadrature one needs an error estimate.  A natural approach is
to evaluate two rules of different orders and compare results.  The
Gauss-Kronrod strategy does this efficiently:

1. Start with an N-point Gauss rule (exact for polynomials of degree $2N-1$).
2. Add N+1 **Kronrod** nodes to form a 2N+1 point rule.
3. The 2N+1 point rule reuses all N Gauss nodes — no extra function evaluations are wasted.

The advantage of using a Gauss and Gauss-Kronrod pair is that the second
rule, which uses $2N+1$ points, actually includes the $N$ points in the
previous Gauss rule.  This means that the function values from that
computation can be reused.  This efficiency comes at the cost of a mild
reduction in the degree of polynomial precision of the Gauss-Kronrod rule.

### Algorithm

The Kronrod nodes are chosen to maximise the degree of polynomial precision
of the combined 2N+1 point rule.  The algorithm proceeds as follows:

1. Compute the Chebyshev coefficients of the orthogonal polynomial whose zeros
   are the Kronrod abscissas.
2. Find approximate values for the abscissas using trigonometric identities.
3. Refine each abscissa using Newton's method.
4. Compute the corresponding weights.

### Storage

Given $N$, let $M = (N+1)/2$.  The Gauss-Kronrod rule will include $2N+1$
points.  By symmetry, only $N+1$ of them need to be stored.

The arrays $x$, $w_1$ and $w_2$ contain the non-negative abscissas in
decreasing order, and the Gauss-Kronrod and Gauss weights respectively.
About half the entries in $w_2$ are zero (corresponding to Kronrod-only
points).

### Example: N = 3

| $i$ | $x$ | $w_1$ (Gauss-Kronrod) | $w_2$ (Gauss) |
|-----|------|----------------------|----------------|
| 1 | 0.960491 | 0.104656 | 0.000000 |
| 2 | 0.774597 | 0.268488 | 0.555556 |
| 3 | 0.434244 | 0.401397 | 0.000000 |
| 4 | 0.000000 | 0.450917 | 0.888889 |

### Example: N = 4

| $i$ | $x$ | $w_1$ (Gauss-Kronrod) | $w_2$ (Gauss) |
|-----|------|----------------------|----------------|
| 1 | 0.976560 | 0.062977 | 0.000000 |
| 2 | 0.861136 | 0.170054 | 0.347855 |
| 3 | 0.640286 | 0.266798 | 0.000000 |
| 4 | 0.339981 | 0.326949 | 0.652145 |
| 5 | 0.000000 | 0.346443 | 0.000000 |

---

## Interval adjustment

The quadrature rule computed by `kronrod` is for the interval $[-1, +1]$.
To integrate over a general interval $[a, b]$, the abscissas and weights
are linearly mapped:

$$
x_i' = \frac{(1 - x_i)\,a + (1 + x_i)\,b}{2}, \qquad
w_i' = \frac{b - a}{2}\, w_i.
$$

The `kronrod_adjust` function performs this transformation.

---

## Error estimation

Given the Gauss estimate $I_G$ (using weights $w_2$) and the Gauss-Kronrod
estimate $I_{GK}$ (using weights $w_1$), the error is estimated as:

$$
\varepsilon \approx |I_{GK} - I_G|.
$$

When this difference is small enough, the integral is considered converged.
