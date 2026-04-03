/**
 * @file kronrod.c
 * @brief Gauss-Kronrod quadrature rule computation.
 *
 * Compute the abscissas and weights of the 2N+1 point Gauss-Kronrod
 * quadrature formula obtained by optimally adding N+1 points to an
 * N-point Gauss quadrature rule.  The advantage of using a Gauss and
 * Gauss-Kronrod pair is that the second rule, which uses 2*N+1 points,
 * actually includes the N points in the previous Gauss rule so that the
 * function values from that computation can be reused.
 *
 * @author  Robert Piessens, Maria Branders (original FORTRAN77)
 * @author  John Burkardt (C translation)
 * @date    2010
 * @note    Distributed under the GNU LGPL license.
 *
 * @par Reference
 *   Robert Piessens, Maria Branders,
 *   "A Note on the Optimal Addition of Abscissas to Quadrature Formulas
 *   of Gauss and Lobatto",
 *   Mathematics of Computation,
 *   Volume 28, Number 125, January 1974, pages 135-139.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "kronrod.h"

/**
 * @brief Calculate a Kronrod abscissa and weight.
 *
 * Uses Newton iteration to refine an initial estimate of a Kronrod
 * abscissa, then computes the corresponding weight.
 *
 * @param n      Order of the Gauss rule.
 * @param m      The value of (N + 1) / 2.
 * @param eps    Requested absolute accuracy of the abscissas.
 * @param coef2  A value needed to compute weights.
 * @param even   Non-zero if N is even.
 * @param b      Chebyshev coefficients, length M+1.
 * @param[in,out] x  On input an estimate; on output the computed abscissa.
 * @param[out]    w  The Kronrod weight.
 */
void abwe1(int n, int m, double eps, double coef2, int even,
           double b[], double *x, double *w)
{
    double ai, b0, b1, b2;
    double d0, d1, d2;
    double delta, dif, f, fd;
    double yy;
    int i, iter, k, ka;

    ka = (*x == 0.0) ? 1 : 0;

    /* Iterative process for the computation of a Kronrod abscissa. */
    for (iter = 1; iter <= 50; iter++) {
        b1 = 0.0;
        b2 = b[m];
        yy = 4.0 * (*x) * (*x) - 2.0;
        d1 = 0.0;

        if (even) {
            ai = m + m + 1;
            d2 = ai * b[m];
            dif = 2.0;
        } else {
            ai = m + 1;
            d2 = 0.0;
            dif = 1.0;
        }

        for (k = 1; k <= m; k++) {
            ai = ai - dif;
            i = m - k + 1;
            b0 = b1;
            b1 = b2;
            d0 = d1;
            d1 = d2;
            b2 = yy * b1 - b0 + b[i - 1];
            if (!even)
                i = i + 1;
            d2 = yy * d1 - d0 + ai * b[i - 1];
        }

        if (even) {
            f = (*x) * (b2 - b1);
            fd = d2 + d1;
        } else {
            f = 0.5 * (b2 - b0);
            fd = 4.0 * (*x) * d2;
        }

        /* Newton correction. */
        delta = f / fd;
        *x = *x - delta;

        if (ka == 1)
            break;

        if (r8_abs(delta) <= eps)
            ka = 1;
    }

    /* Catch non-convergence. */
    if (ka != 1) {
        fprintf(stderr, "\nabwe1 - Fatal error!\n");
        fprintf(stderr, "  Iteration limit reached.\n");
        fprintf(stderr, "  eps was %e\n", eps);
        fprintf(stderr, "  Last delta was %e\n", delta);
        exit(1);
    }

    /* Computation of the weight. */
    d0 = 1.0;
    d1 = *x;
    ai = 0.0;
    for (k = 2; k <= n; k++) {
        ai = ai + 1.0;
        d2 = ((ai + ai + 1.0) * (*x) * d1 - ai * d0) / (ai + 1.0);
        d0 = d1;
        d1 = d2;
    }

    *w = coef2 / (fd * d2);
}

/**
 * @brief Calculate a Gaussian abscissa and two weights.
 *
 * Uses Newton iteration to refine an initial estimate of a Gaussian
 * abscissa, then computes the Gauss-Kronrod and Gauss weights.
 *
 * @param n      Order of the Gauss rule.
 * @param m      The value of (N + 1) / 2.
 * @param eps    Requested absolute accuracy of the abscissas.
 * @param coef2  A value needed to compute weights.
 * @param even   Non-zero if N is even.
 * @param b      Chebyshev coefficients, length M+1.
 * @param[in,out] x   On input an estimate; on output the computed abscissa.
 * @param[out]    w1   The Gauss-Kronrod weight.
 * @param[out]    w2   The Gauss weight.
 */
void abwe2(int n, int m, double eps, double coef2, int even,
           double b[], double *x, double *w1, double *w2)
{
    double ai, an, delta;
    double p0, p1, p2;
    double pd0, pd1, pd2;
    double yy;
    int i, iter, k, ka;

    ka = (*x == 0.0) ? 1 : 0;

    /* Iterative process for the computation of a Gaussian abscissa. */
    for (iter = 1; iter <= 50; iter++) {
        p0 = 1.0;
        p1 = *x;
        pd0 = 0.0;
        pd1 = 1.0;

        /*
         * When N is 1, we need to initialise P2 and PD2 to avoid
         * problems with DELTA.
         */
        if (n <= 1) {
            if (r8_epsilon() < r8_abs(*x)) {
                p2 = (3.0 * (*x) * (*x) - 1.0) / 2.0;
                pd2 = 3.0 * (*x);
            } else {
                p2 = 3.0 * (*x);
                pd2 = 3.0;
            }
        }

        ai = 0.0;
        for (k = 2; k <= n; k++) {
            ai = ai + 1.0;
            p2 = ((ai + ai + 1.0) * (*x) * p1 - ai * p0) / (ai + 1.0);
            pd2 = ((ai + ai + 1.0) * (p1 + (*x) * pd1) - ai * pd0)
                  / (ai + 1.0);
            p0 = p1;
            p1 = p2;
            pd0 = pd1;
            pd1 = pd2;
        }

        /* Newton correction. */
        delta = p2 / pd2;
        *x = *x - delta;

        if (ka == 1)
            break;

        if (r8_abs(delta) <= eps)
            ka = 1;
    }

    /* Catch non-convergence. */
    if (ka != 1) {
        fprintf(stderr, "\nabwe2 - Fatal error!\n");
        fprintf(stderr, "  Iteration limit reached.\n");
        fprintf(stderr, "  eps was %e\n", eps);
        fprintf(stderr, "  Last delta was %e\n", delta);
        exit(1);
    }

    /* Computation of the weight. */
    an = n;
    *w2 = 2.0 / (an * pd2 * p0);

    p1 = 0.0;
    p2 = b[m];
    yy = 4.0 * (*x) * (*x) - 2.0;
    for (k = 1; k <= m; k++) {
        i = m - k + 1;
        p0 = p1;
        p1 = p2;
        p2 = yy * p1 - p0 + b[i - 1];
    }

    if (even)
        *w1 = *w2 + coef2 / (pd2 * (*x) * (p2 - p1));
    else
        *w1 = *w2 + 2.0 * coef2 / (pd2 * (p2 - p0));
}

/**
 * @brief Add N+1 Kronrod points to an N-point Gaussian rule.
 *
 * Compute the abscissas and weights of the 2N+1 point Gauss-Kronrod
 * quadrature formula for integration over [-1, +1].  Only the
 * non-negative abscissas are returned (the rule is symmetric).
 *
 * Given N, let M = (N + 1) / 2.  The arrays @p x, @p w1 and @p w2
 * contain the non-negative abscissas in decreasing order, and the
 * weights of each abscissa in the Gauss-Kronrod and Gauss rules
 * respectively.  About half the entries in @p w2 are zero.
 *
 * @param n    Order of the Gauss rule.
 * @param eps  Requested absolute accuracy of the abscissas.
 * @param[out] x   Abscissas, length N+1 (non-negative, descending).
 * @param[out] w1  Gauss-Kronrod weights, length N+1.
 * @param[out] w2  Gauss weights, length N+1 (zero for Kronrod-only points).
 */
void kronrod(int n, double eps, double x[], double w1[], double w2[])
{
    double ak, an, bb, c, coef, coef2, d, s;
    double x1, xx, y;
    double *b, *tau;
    int even, i, k, l, ll, m;

    b = (double *)malloc((size_t)(((n + 1) / 2) + 1) * sizeof(double));
    tau = (double *)malloc((size_t)((n + 1) / 2) * sizeof(double));

    m = (n + 1) / 2;
    even = (2 * m == n);

    d = 2.0;
    an = 0.0;
    for (k = 1; k <= n; k++) {
        an = an + 1.0;
        d = d * an / (an + 0.5);
    }

    /* Calculation of the Chebyshev coefficients of the orthogonal polynomial. */
    tau[0] = (an + 2.0) / (an + an + 3.0);
    b[m - 1] = tau[0] - 1.0;
    ak = an;

    for (l = 1; l < m; l++) {
        ak = ak + 2.0;
        tau[l] = ((ak - 1.0) * ak
                  - an * (an + 1.0)) * (ak + 2.0) * tau[l - 1]
                 / (ak * ((ak + 3.0) * (ak + 2.0)
                           - an * (an + 1.0)));
        b[m - l - 1] = tau[l];

        for (ll = 1; ll <= l; ll++)
            b[m - l - 1] = b[m - l - 1] + tau[ll - 1] * b[m - l + ll - 1];
    }

    b[m] = 1.0;

    /* Calculation of approximate values for the abscissas. */
    bb = sin(1.570796 / (an + an + 1.0));
    x1 = sqrt(1.0 - bb * bb);
    s = 2.0 * bb * x1;
    c = sqrt(1.0 - s * s);
    coef = 1.0 - (1.0 - 1.0 / an) / (8.0 * an * an);
    xx = coef * x1;

    /*
     * Coefficient needed for weights.
     * COEF2 = 2^(2*n+1) * n! * n! / (2n+1)!
     */
    coef2 = 2.0 / (double)(2 * n + 1);
    for (i = 1; i <= n; i++)
        coef2 = coef2 * 4.0 * (double)(i) / (double)(n + i);

    /*
     * Calculation of the K-th abscissa (a Kronrod abscissa)
     * and the corresponding weight.
     */
    for (k = 1; k <= n; k = k + 2) {
        abwe1(n, m, eps, coef2, even, b, &xx, w1 + k - 1);
        w2[k - 1] = 0.0;

        x[k - 1] = xx;
        y = x1;
        x1 = y * c - bb * s;
        bb = y * s + bb * c;

        if (k == n)
            xx = 0.0;
        else
            xx = coef * x1;

        /*
         * Calculation of the K+1 abscissa (a Gaussian abscissa)
         * and the corresponding weights.
         */
        abwe2(n, m, eps, coef2, even, b, &xx, w1 + k, w2 + k);

        x[k] = xx;
        y = x1;
        x1 = y * c - bb * s;
        bb = y * s + bb * c;
        xx = coef * x1;
    }

    /*
     * If N is even, we have one more Kronrod abscissa to compute,
     * namely the origin.
     */
    if (even) {
        xx = 0.0;
        abwe1(n, m, eps, coef2, even, b, &xx, w1 + n);
        w2[n] = 0.0;
        x[n] = xx;
    }

    free(b);
    free(tau);
}

/**
 * @brief Adjust a Gauss-Kronrod rule from [-1, +1] to [a, b].
 *
 * Linearly maps the abscissas and scales the weights so that the
 * rule integrates over the interval [a, b] instead of [-1, +1].
 *
 * @param a   Left endpoint of the new interval.
 * @param b   Right endpoint of the new interval.
 * @param n   Order of the Gauss rule.
 * @param[in,out] x   Abscissas, length N+1.
 * @param[in,out] w1  Gauss-Kronrod weights, length N+1.
 * @param[in,out] w2  Gauss weights, length N+1.
 */
void kronrod_adjust(double a, double b, int n,
                    double x[], double w1[], double w2[])
{
    int i;

    for (i = 0; i < n + 1; i++) {
        x[i]  = ((1.0 - x[i]) * a + (1.0 + x[i]) * b) / 2.0;
        w1[i] = ((b - a) / 2.0) * w1[i];
        w2[i] = ((b - a) / 2.0) * w2[i];
    }
}

/**
 * @brief Return the absolute value of a double.
 *
 * @param x  Input value.
 * @return   Absolute value of @p x.
 */
double r8_abs(double x)
{
    return (0.0 <= x) ? x : -x;
}

/**
 * @brief Return the machine epsilon for double precision.
 *
 * R8_EPSILON is a number R which is a power of 2 with the property
 * that, to the precision of the computer's arithmetic, 1 < 1 + R
 * but 1 = (1 + R / 2).
 *
 * @return  The smallest double e such that 1 + e > 1.
 */
double r8_epsilon(void)
{
    static double value = 2.220446049250313E-016;

    return value;
}

/**
 * @brief Print the current date and time as a timestamp.
 *
 * Example output: "31 May 2001 09:45:54 AM".
 */
void timestamp(void)
{
#define TIME_SIZE 40
    static char time_buffer[TIME_SIZE];
    const struct tm *tm_ptr;
    time_t now;

    now = time(NULL);
    tm_ptr = localtime(&now);
    strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr);
    printf("%s\n", time_buffer);
#undef TIME_SIZE
}
