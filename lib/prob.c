#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "bug.h"
#include "prob.h"

int prob_cdf_init(struct prob_cdf *cdf, const size_t len)
{
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
 * prob_cdf_inverse - sample from the CDF via inverse
 * transform.
 *
 * Hybrid binary/linear search.  The binary phase repeatedly
 * halves the stride d and advances i only when F(i+d) <= y,
 * converging to the neighbourhood of the target in O(log n).
 * The linear scan handles the remaining entries where the
 * binary stride has reached zero.  Total cost: O(log n + k)
 * where k is a small constant from the linear tail.
 */
size_t prob_cdf_inverse(const struct prob_cdf *cdf, const double y)
{
	size_t i = 0, d = cdf->len;

	while ((d /= 2) > 0) {
		BUG_ON(i >= cdf->len);
		if (cdf->y[i + d] <= y)
			i += d;
	}
	while (i < cdf->len - 1 && cdf->y[i] <= y)
		i++;

	return i;
}
