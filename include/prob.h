#ifndef PROB_H
#define PROB_H

/*
 * Discrete CDF builder and inverse-sampler over
 * possibly-negative weighted samples.
 *
 * Two operations:
 *   - build:  normalise |w[i]| over n samples into
 *             a cumulative distribution F(i).
 *   - lookup: given y in [0, 1], return the index
 *             sampled with probability f(i) =
 *             F(i) - F(i-1).  Sampler convention
 *             (smallest i with F(i) > y, clamped).
 *
 * No thread-safety guarantees: each struct prob_cdf
 * is owned by one writer at a time; concurrent
 * readers of an unchanging CDF are fine.  Used by
 * the qDRIFT-family samplers in circ/.
 */

#include <stddef.h>

struct prob_cdf {
	double *y;
	size_t len;
};

/* Allocate cdf->y for `len` entries and record the
 * length.  Returns -1 if `len == 0` or the malloc
 * fails (the `len == 0` rejection is an enforced
 * contract -- the rest of the API assumes len > 0),
 * 0 on success.  The buffer is uninitialised; call
 * prob_cdf_from_array_strided to fill it. */
int prob_cdf_init(struct prob_cdf *cdf, size_t len);

/* Release the buffer and zero cdf->y + cdf->len so
 * a second free is a clean no-op.  Safe to call on
 * a default-zero struct prob_cdf. */
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
