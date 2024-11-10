#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "bug.h"
#include "phase2/prob.h"

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

int prob_cdf_from(struct prob_cdf *cdf, double (*get_vals)(void *), void *data)
{
	/* Calculate PDF. */
	double lambda = 0.0;
	for (double *y = cdf->y; y < cdf->y + cdf->len; y++) {
		const double yi = fabs(get_vals(data));
		*y = yi;
		lambda += yi;
	}
	if (lambda < DBL_EPSILON)
		return -1;

	/* Calculate CDF */
	double f = 0.0;
	for (double *y = cdf->y; y < cdf->y + cdf->len; y++) {
		f += *y / lambda;
		*y = f;
	}

	return 0;
}

size_t prob_cdf_inverse(const struct prob_cdf *cdf, const double y)
{
	size_t i = 0, d = cdf->len;

	while ((d /= 2) > 0) {
		BUG_ON(i >= cdf->len);
		if (cdf->y[i + d] <= y)
			i += d;
	}
	while (i < cdf->len && cdf->y[i] <= y)
		i++;

	return i;
}
