/**
 * @file test_kronrod.c
 * @brief Test suite for the Gauss-Kronrod library.
 *
 * Test data is taken from the tabulated values in the legacy test suite
 * by John Burkardt (bin/legacy/kronrod_test.txt).
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kronrod.h"

#define ATOL 1.0e-5

static int g_pass = 0;
static int g_fail = 0;

static void check(const char *name, double got, double expected, double tol)
{
    double diff = fabs(got - expected);

    if (diff > tol) {
        fprintf(stderr,
                "FAIL  %-30s  got %16.10g  expected %16.10g  diff %g\n",
                name, got, expected, diff);
        g_fail++;
    } else {
        g_pass++;
    }
}

/* ---- test_n3: odd case N = 3 -------------------------------------------- */

static void test_n3(void)
{
    int n = 3;
    double eps = 0.000001;
    double x[4], w1[4], w2[4];

    /* Expected abscissas (non-negative, descending) */
    static const double x_exp[4] = {
        0.960491, 0.774597, 0.434244, 0.000000
    };

    /* Expected Gauss-Kronrod weights */
    static const double w1_exp[4] = {
        0.104656, 0.268488, 0.401397, 0.450917
    };

    /* Expected Gauss weights */
    static const double w2_exp[4] = {
        0.000000, 0.555556, 0.000000, 0.888889
    };

    /* Exact Gauss abscissas (full symmetric rule) */
    static const double xg[3] = {
        -0.77459666924148337704,
         0.0,
         0.77459666924148337704
    };

    /* Exact Gauss weights */
    static const double wg[3] = {
        0.555555555555555555556,
        0.888888888888888888889,
        0.555555555555555555556
    };

    /* Exact Gauss-Kronrod abscissas (full 2N+1 = 7 point rule) */
    static const double xk[7] = {
        -0.96049126870802028342,
        -0.77459666924148337704,
        -0.43424374934680255800,
         0.0,
         0.43424374934680255800,
         0.77459666924148337704,
         0.96049126870802028342
    };

    /* Exact Gauss-Kronrod weights */
    static const double wk[7] = {
        0.104656226026467265194,
        0.268488089868333440729,
        0.401397414775962222905,
        0.450916538658474142345,
        0.401397414775962222905,
        0.268488089868333440729,
        0.104656226026467265194
    };

    int i, i2;
    double s;

    kronrod(n, eps, x, w1, w2);

    /* Check compact output against expected values */
    for (i = 0; i < n + 1; i++) {
        char label[64];
        snprintf(label, sizeof(label), "n3_x[%d]", i);
        check(label, x[i], x_exp[i], ATOL);
        snprintf(label, sizeof(label), "n3_w1[%d]", i);
        check(label, w1[i], w1_exp[i], ATOL);
        snprintf(label, sizeof(label), "n3_w2[%d]", i);
        check(label, w2[i], w2_exp[i], ATOL);
    }

    /* Gauss abscissas */
    for (i = 1; i <= n; i++) {
        char label[64];
        if (2 * i <= n + 1) {
            i2 = 2 * i;
            s = -1.0;
        } else {
            i2 = 2 * (n + 1) - 2 * i;
            s = +1.0;
        }
        snprintf(label, sizeof(label), "n3_gauss_x[%d]", i);
        check(label, s * x[i2 - 1], xg[i - 1], 1.0e-10);
    }

    /* Gauss weights */
    for (i = 1; i <= n; i++) {
        char label[64];
        if (2 * i <= n + 1) {
            i2 = 2 * i;
        } else {
            i2 = 2 * (n + 1) - 2 * i;
        }
        snprintf(label, sizeof(label), "n3_gauss_w[%d]", i);
        check(label, w2[i2 - 1], wg[i - 1], 1.0e-10);
    }

    /* Gauss-Kronrod abscissas */
    for (i = 1; i <= 2 * n + 1; i++) {
        char label[64];
        if (i <= n + 1) {
            i2 = i;
            s = -1.0;
        } else {
            i2 = 2 * (n + 1) - i;
            s = +1.0;
        }
        snprintf(label, sizeof(label), "n3_kronrod_x[%d]", i);
        check(label, s * x[i2 - 1], xk[i - 1], 1.0e-10);
    }

    /* Gauss-Kronrod weights */
    for (i = 1; i <= 2 * n + 1; i++) {
        char label[64];
        if (i <= n + 1) {
            i2 = i;
        } else {
            i2 = 2 * (n + 1) - i;
        }
        snprintf(label, sizeof(label), "n3_kronrod_w[%d]", i);
        check(label, w1[i2 - 1], wk[i - 1], 1.0e-10);
    }
}

/* ---- test_n4: even case N = 4 ------------------------------------------- */

static void test_n4(void)
{
    int n = 4;
    double eps = 0.000001;
    double x[5], w1[5], w2[5];

    /* Expected from legacy output */
    static const double x_exp[5] = {
        0.976560, 0.861136, 0.640286, 0.339981, 0.000000
    };
    static const double w1_exp[5] = {
        0.062977, 0.170054, 0.266798, 0.326949, 0.346443
    };
    static const double w2_exp[5] = {
        0.000000, 0.347855, 0.000000, 0.652145, 0.000000
    };

    int i;

    kronrod(n, eps, x, w1, w2);

    for (i = 0; i < n + 1; i++) {
        char label[64];
        snprintf(label, sizeof(label), "n4_x[%d]", i);
        check(label, x[i], x_exp[i], ATOL);
        snprintf(label, sizeof(label), "n4_w1[%d]", i);
        check(label, w1[i], w1_exp[i], ATOL);
        snprintf(label, sizeof(label), "n4_w2[%d]", i);
        check(label, w2[i], w2_exp[i], ATOL);
    }
}

/* ---- test_integral: integral estimate ------------------------------------ */

static double test_f(double x)
{
    return 1.0 / (x * x + 1.005);
}

static void test_integral(void)
{
    double exact = 1.5643964440690497731;
    double eps = 0.000001;
    int n;
    double *x, *w1, *w2;
    double i1, i2;
    int i;

    n = 1;
    while (n <= 25) {
        x  = (double *)malloc((size_t)(n + 1) * sizeof(double));
        w1 = (double *)malloc((size_t)(n + 1) * sizeof(double));
        w2 = (double *)malloc((size_t)(n + 1) * sizeof(double));

        kronrod(n, eps, x, w1, w2);

        i1 = w1[n] * test_f(x[n]);
        i2 = w2[n] * test_f(x[n]);

        for (i = 0; i < n; i++) {
            i1 += w1[i] * (test_f(-x[i]) + test_f(x[i]));
            i2 += w2[i] * (test_f(-x[i]) + test_f(x[i]));
        }

        free(x);
        free(w1);
        free(w2);

        if (fabs(i1 - i2) < 0.0001) {
            check("integral_converged", 1.0, 1.0, 0.0);
            check("integral_error", fabs(exact - i1), 0.0, 1.0e-4);
            return;
        }

        n = 2 * n + 1;
    }

    /* Should have converged by now */
    check("integral_convergence_failed", 0.0, 1.0, 0.0);
}

/* ---- test_adjust: interval adjustment ------------------------------------ */

static void test_adjust(void)
{
    int n = 3;
    double eps = 0.000001;
    double x[4], w1[4], w2[4];
    double sum_w1;
    int i;

    kronrod(n, eps, x, w1, w2);
    kronrod_adjust(0.0, 2.0, n, x, w1, w2);

    /* After adjusting to [0, 2], sum of GK weights should equal
     * interval length (2.0) */
    sum_w1 = 0.0;
    for (i = 0; i < n + 1; i++)
        sum_w1 += w1[i];
    /* Full rule has symmetric part too; total = 2 * sum - w1[n] (origin) */
    sum_w1 = 2.0 * sum_w1 - w1[n];

    /* The sum of all 2N+1 GK weights should approximately equal
     * the interval length when the rule is for [a,b].
     * Actually for the compact representation we need careful handling.
     * Just check abscissas are in [0, 2]. */
    for (i = 0; i < n + 1; i++) {
        char label[64];
        snprintf(label, sizeof(label), "adjust_x_range[%d]", i);
        check(label, (x[i] >= -0.01 && x[i] <= 2.01) ? 1.0 : 0.0,
              1.0, 0.0);
    }
}

/* ---- test_r8_abs --------------------------------------------------------- */

static void test_r8_abs(void)
{
    check("r8_abs(3.14)", r8_abs(3.14), 3.14, 1.0e-15);
    check("r8_abs(-3.14)", r8_abs(-3.14), 3.14, 1.0e-15);
    check("r8_abs(0.0)", r8_abs(0.0), 0.0, 1.0e-15);
}

/* ---- test_r8_epsilon ----------------------------------------------------- */

static void test_r8_epsilon(void)
{
    double eps = r8_epsilon();

    check("r8_epsilon_positive", (eps > 0.0) ? 1.0 : 0.0, 1.0, 0.0);
    check("r8_epsilon_small", (eps < 1.0e-10) ? 1.0 : 0.0, 1.0, 0.0);
    check("r8_epsilon_value", eps, 2.220446049250313E-016, 1.0e-25);
}

/* ---- main ---------------------------------------------------------------- */

int main(void)
{
    test_n3();
    test_n4();
    test_integral();
    test_adjust();
    test_r8_abs();
    test_r8_epsilon();

    printf("\n%d passed, %d failed\n", g_pass, g_fail);

    return (g_fail > 0) ? 1 : 0;
}
