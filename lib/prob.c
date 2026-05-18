/*
 * lib/prob.c -- discrete CDF builder and
 * inverse-sampler.  See include/prob.h for the public
 * API and contract.
 *
 * The module is used by the qDRIFT (circ/qdrift.c)
 * and composite (circ/cmpsit.c) samplers; both build
 * a CDF over the Hamiltonian's |c_k| weights at init
 * time and call prob_cdf_inverse(cdf, x) per random
 * draw to map a uniform x in [0, 1] to a term index.
 *
 * Implementation notes:
 *   - prob_cdf_from_array_strided is a two-pass
 *     normalisation; lambda (the L1 norm) is exposed
 *     to callers via the `out_lambda` parameter.
 *   - prob_cdf_inverse uses a hybrid binary/linear
 *     walk; bench/b-prob compares it against a
 *     textbook upper_bound and the hybrid wins by
 *     ~20% at the typical n in [100, 1000].
 *
 * No dependencies beyond libc.
 */

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "prob.h"

int prob_cdf_init(struct prob_cdf *cdf, const size_t len)
{
	if (len == 0)
		return -1;

	double *x = malloc(sizeof *x * len);
	if (!x)
		return -1;

	cdf->y = x;
	cdf->len = len;

	return 0;
}

void prob_cdf_free(struct prob_cdf *cdf)
{
	free(cdf->y);
	cdf->y = NULL;
	cdf->len = 0;
}

/*
 * prob_cdf_from_array_strided - build a CDF from
 * cdf->len possibly-negative weights at `stride`-byte
 * intervals from `base`.
 *
 * Two-pass construction:
 *  1. Accumulate |w[i]| into cdf->y[i] and sum into lambda.
 *  2. Normalise each entry by lambda and compute the
 *     running sum to form the CDF: y[i] = sum_{k<=i} p(k).
 *
 * Returns -1 if lambda < DBL_EPSILON (all weights
 * negligible).  Writes lambda through `out_lambda` on
 * success when non-NULL.
 */
int prob_cdf_from_array_strided(struct prob_cdf *cdf,
	const double *base, size_t stride, double *out_lambda)
{
	double lambda = 0.0;
	for (size_t i = 0; i < cdf->len; i++) {
		const double *w = (const double *)
			((const char *)base + i * stride);
		const double yi = fabs(*w);
		cdf->y[i] = yi;
		lambda += yi;
	}
	if (lambda < DBL_EPSILON)
		return -1;

	double f = 0.0;
	for (size_t i = 0; i < cdf->len; i++) {
		f += cdf->y[i] / lambda;
		cdf->y[i] = f;
	}

	if (out_lambda)
		*out_lambda = lambda;
	return 0;
}

/*
 * prob_cdf_inverse - inverse-CDF lookup (sampler
 * convention).
 *
 * Returns the smallest index i such that cdf->y[i] > y,
 * or cdf->len - 1 when no such index exists.  For a
 * uniform draw y in [0, 1] this samples index i with
 * probability f(i) = F(i) - F(i-1).
 *
 * Hybrid binary/linear walk: the binary phase halves
 * `d` and advances `i` only when y[i+d] <= y, converging
 * to the neighbourhood of the answer in O(log n).  The
 * linear tail mops up the few entries left over when
 * cdf->len is not a power of two.  Faster than a
 * textbook upper_bound binary search at typical sizes
 * (n in the hundreds) because each iteration carries
 * less arithmetic.
 */
size_t prob_cdf_inverse(const struct prob_cdf *cdf, const double y)
{
	size_t i = 0, d = cdf->len;

	while ((d /= 2) > 0)
		if (cdf->y[i + d] <= y)
			i += d;
	while (i < cdf->len - 1 && cdf->y[i] <= y)
		i++;

	return i;
}
