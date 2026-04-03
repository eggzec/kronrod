/**
 * @file kronrod.h
 * @brief Gauss-Kronrod quadrature rule computation.
 *
 * Compute the abscissas and weights of the 2N+1 point Gauss-Kronrod
 * quadrature formula obtained by optimally adding N+1 points to an
 * N-point Gauss quadrature rule.
 *
 * @author  Robert Piessens, Maria Branders (original FORTRAN77)
 * @author  John Burkardt (C translation)
 * @date    2010
 *
 * @par License
 *   Distributed under the GNU LGPL license.
 *
 * @par Reference
 *   Robert Piessens, Maria Branders,
 *   "A Note on the Optimal Addition of Abscissas to Quadrature Formulas
 *   of Gauss and Lobatto",
 *   Mathematics of Computation,
 *   Volume 28, Number 125, January 1974, pages 135-139.
 */

#ifndef KRONROD_H
#define KRONROD_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Calculate a Kronrod abscissa and weight.
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
           double b[], double *x, double *w);

/**
 * @brief Calculate a Gaussian abscissa and two weights.
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
           double b[], double *x, double *w1, double *w2);

/**
 * @brief Add N+1 Kronrod points to an N-point Gaussian rule.
 *
 * Compute the abscissas and weights of the 2N+1 point Gauss-Kronrod
 * quadrature formula for integration over [-1, +1].  Only the
 * non-negative abscissas are returned (the rule is symmetric).
 *
 * @param n    Order of the Gauss rule.
 * @param eps  Requested absolute accuracy of the abscissas.
 * @param[out] x   Abscissas, length N+1 (non-negative, descending).
 * @param[out] w1  Gauss-Kronrod weights, length N+1.
 * @param[out] w2  Gauss weights, length N+1 (zero for Kronrod-only points).
 */
void kronrod(int n, double eps, double x[], double w1[], double w2[]);

/**
 * @brief Adjust a Gauss-Kronrod rule from [-1, +1] to [a, b].
 *
 * @param a   Left endpoint of the new interval.
 * @param b   Right endpoint of the new interval.
 * @param n   Order of the Gauss rule.
 * @param[in,out] x   Abscissas, length N+1.
 * @param[in,out] w1  Gauss-Kronrod weights, length N+1.
 * @param[in,out] w2  Gauss weights, length N+1.
 */
void kronrod_adjust(double a, double b, int n,
                    double x[], double w1[], double w2[]);

/**
 * @brief Return the absolute value of a double.
 *
 * @param x  Input value.
 * @return   Absolute value of @p x.
 */
double r8_abs(double x);

/**
 * @brief Return the machine epsilon for double precision.
 *
 * @return  The smallest double e such that 1 + e > 1.
 */
double r8_epsilon(void);

/**
 * @brief Print the current date and time as a timestamp.
 */
void timestamp(void);

#ifdef __cplusplus
}
#endif

#endif /* KRONROD_H */
