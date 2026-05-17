#ifndef PROB_H
#define PROB_H

#include <stddef.h>

struct prob_cdf {
	double *y;
	size_t len;
};

/* Allocate the CDF buffer.  `len` must be greater than 0;
 * returns -1 if it is 0 or the allocation fails, else 0. */
int prob_cdf_init(struct prob_cdf *cdf, size_t len);

/* Release the buffer and zero `y` + `len` so a second
 * free is a clean no-op. */
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
 * Inverse-CDF lookup (sampler convention).
 *
 * Returns the smallest index i such that F(i) > y, or
 * cdf->len - 1 when no such index exists.  For a uniform
 * draw y in [0, 1] this samples index i with probability
 * f(i) = F(i) - F(i-1).
 */
size_t prob_cdf_inverse(const struct prob_cdf *cdf, double y);

#endif // PROB_H
