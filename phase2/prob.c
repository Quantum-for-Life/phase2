#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "phase2/prob.h"

int prob_pd_init(struct prob_pd *pd, const size_t len)
{
	double *x = malloc(sizeof *x * len);
	if (!x)
		return -1;

	pd->x = x;
	pd->len = len;

	return 0;
}

void prob_pd_free(struct prob_pd *pd)
{
	free(pd->x);
}

int prob_pdf_from_samples(
	struct prob_pd *pd, double (*get_smpl)(void *), void *data)
{
	double lambda = 0.0;
	for (size_t i = 0; i < pd->len; i++) {
		const double xi = fabs(get_smpl(data));
		pd->x[i] = xi;
		lambda += xi;
	}
	if (lambda < DBL_EPSILON)
		return -1;

	for (double *x = pd->x; x < pd->x + pd->len; x++)
		*x /= lambda;

	return 0;
}

size_t prob_cdf_inverse(const struct prob_pd *pd, const double y)
{
	size_t i = 0;
	double cdf = 0.0;

	while (cdf <= y)
		cdf += pd->x[i++];

	return i - 1;
}
