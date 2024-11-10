#ifndef PROB_H
#define PROB_H

#include <stddef.h>

struct prob_pd {
	double *x; /* normalized array of nonnegative doubles */
	size_t len;
};

int prob_pd_init(struct prob_pd *pd, size_t len);

void prob_pd_free(struct prob_pd *pd);

/*
 * Calculate probability distribution (probability density function) from
 * a set of not normalized, possibly negative samples of length pd->len.
 *
 * The pdf is calculated by calling get_smpl() pd->len times, taking
 * the absolute value, and finally normalizing.
 */
int prob_pdf_from_samples(
	struct prob_pd *pd, double (*get_smpl)(void *), void *data);

/*
 * Compute the inverse of a discrete cumulative distribution function (CDF).
 *
 * The CDF is defined as:
 *     F(i) = \sum_{j < i} f(j),
 *   where f(j) is the probability density function.
 *
 * Arguments:
 *   pdf - probability density function (consisting of a finite set of
 *         nonnegative numbers)
 *   y   - argument of F^{-1}(y)
 *
 * Returns:
 *       Index i of the corresponding bin such that F(i) = y.
 */
size_t prob_cdf_inverse(const struct prob_pd *pd, double y);

#endif // PROB_H
