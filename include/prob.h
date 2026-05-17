#ifndef PROB_H
#define PROB_H

#include <stddef.h>

struct prob_cdf {
	double *y;
	size_t len;
};

/* Len must be greater than 0. */
int prob_cdf_init(struct prob_cdf *cdf, size_t len);

void prob_cdf_free(struct prob_cdf *cdf);

/*
 * Build a CDF from cdf->len possibly-negative weights laid out at
 * stride-byte intervals starting at `base` (so weight i is the double
 * at byte offset i*stride).  The absolute values are normalised into
 * a PDF and accumulated into cdf->y.
 *
 * If `out_lambda` is non-NULL, the sum of |w[i]| (the L1 norm of the
 * raw weights) is written there.  Pass NULL to discard.
 *
 * Returns -1 if the total mass is below DBL_EPSILON (the distribution
 * is effectively zero), 0 on success.
 */
int prob_cdf_from_array_strided(struct prob_cdf *cdf,
	const double *base, size_t stride, double *out_lambda);

/*
 * Compute the inverse of a discrete cumulative distribution function (CDF).
 *
 * The CDF is defined as:
 *     F(i) = \sum_{j <= i} f(j),
 *   where f(j) is the probability density function.
 *
 * Arguments:
 *   cdf - pointer to CDF
 *   y   - argument of F^{-1}(y)
 *
 * Returns:
 *       Largest index i such that F(i) <= y, or 0 if F(i) > y for all i.
 */
size_t prob_cdf_inverse(const struct prob_cdf *cdf, double y);

#endif // PROB_H
