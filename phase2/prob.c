#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "phase2/prob.h"

int prob_cdf_init(struct prob_cdf *cdf, const size_t len)
{
	double *x = malloc(sizeof *x * len);
	if (!x)
		return -1;

	cdf->x = x;
	cdf->len = len;

	return 0;
}

void prob_cdf_free(struct prob_cdf *cdf)
{
	free(cdf->x);
}

int prob_cdf_from_samples(
	struct prob_cdf *cdf, double (*get_smpl)(void *), void *data)
{
	double lambda = 0.0;
	for (size_t i = 0; i < cdf->len; i++) {
		const double xi = fabs(get_smpl(data));
		cdf->x[i] = xi;
		lambda += xi;
	}
	if (lambda < DBL_EPSILON)
		return -1;

	for (double *x = cdf->x; x < cdf->x + cdf->len; x++)
		*x /= lambda;

	return 0;
}

size_t prob_cdf_inverse(const struct prob_cdf *cdf, const double y)
{
	size_t i = 0;
	double f = 0.0;

	while (f <= y)
		f += cdf->x[i++];

	return i - 1;
}
