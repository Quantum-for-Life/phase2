#ifndef PROB_H
#define PROB_H

/*
 * Discrete CDF over possibly-negative weights, with
 * inverse-CDF sampling (sampler convention: smallest
 * i with F(i) > y, clamped to len-1).
 */

#include <stddef.h>

struct prob_cdf {
	double *y;
	size_t len;
};

/* Allocate cdf->y for `len > 0` entries.  Returns -1
 * on len == 0 or malloc failure. */
int prob_cdf_init(struct prob_cdf *cdf, size_t len);

/* Free cdf->y; zero cdf->y and cdf->len so repeat
 * free and default-zero structs are safe. */
void prob_cdf_free(struct prob_cdf *cdf);

/* Fill cdf with the CDF of |w[i]| for cdf->len weights
 * at byte offsets i*stride from base.  Writes the L1
 * norm to out_lambda if non-NULL.  Returns -1 if total
 * mass < DBL_EPSILON. */
int prob_cdf_from_array_strided(struct prob_cdf *cdf,
	const double *base, size_t stride, double *out_lambda);

/* Smallest i with cdf->y[i] > y, or cdf->len - 1. */
size_t prob_cdf_inverse(const struct prob_cdf *cdf, double y);

#endif // PROB_H
