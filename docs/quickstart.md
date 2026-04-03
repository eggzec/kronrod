# Quickstart

This page walks through progressively richer examples of using `kronrod`.
Each example is self-contained and can be pasted directly into a Python script
or interactive session. For mathematical background, see the [Theory](theory.md) section.

---

## Basic quadrature rule

Compute a 3-point Gauss / 7-point Gauss-Kronrod rule on [-1, +1]:

```python
import kronrod

x, w1, w2 = kronrod.kronrod(3)

print("Abscissas:", x)
print("Gauss-Kronrod weights:", w1)
print("Gauss weights:", w2)
```

## Even-order rule

Compute a 4-point Gauss / 9-point Gauss-Kronrod rule:

```python
import kronrod

x, w1, w2 = kronrod.kronrod(4)

for i in range(len(x)):
    print(f"  x={x[i]:.6f}  w1={w1[i]:.6f}  w2={w2[i]:.6f}")
```

## Adjusting to a custom interval

Map the rule from [-1, +1] to an arbitrary interval [a, b]:

```python
import kronrod

x, w1, w2 = kronrod.kronrod(3)
x_adj, w1_adj, w2_adj = kronrod.kronrod_adjust(0.0, 1.0, 3, x, w1, w2)

print("Adjusted abscissas:", x_adj)
print("Adjusted Gauss-Kronrod weights:", w1_adj)
```

## Estimating an integral

Use the Gauss-Kronrod pair to estimate an integral with error control:

```python
import math
import kronrod


def f(x):
    return 1.0 / (x * x + 1.005)


exact = 1.5643964440690497731
n = 3

while n <= 25:
    x, w1, w2 = kronrod.kronrod(n)

    # Gauss-Kronrod estimate
    i1 = w1[-1] * f(x[-1])
    for i in range(n):
        i1 += w1[i] * (f(-x[i]) + f(x[i]))

    # Gauss estimate
    i2 = w2[-1] * f(x[-1])
    for i in range(n):
        i2 += w2[i] * (f(-x[i]) + f(x[i]))

    error = abs(i1 - i2)
    print(f"n={n:2d}  GK={i1:.10f}  G={i2:.10f}  err={error:.2e}")

    if error < 0.0001:
        break
    n = 2 * n + 1

print(f"Actual error: {abs(exact - i1):.2e}")
```

## Comparing with SciPy

```python
import kronrod

# kronrod gives the same nodes/weights used internally by
# scipy.integrate.quad for Gauss-Kronrod adaptive quadrature.
x, w1, w2 = kronrod.kronrod(7)
print("7-point Gauss / 15-point Gauss-Kronrod abscissas:")
print(x)
```
