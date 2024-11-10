#ifndef PROB_H
#define PROB_H

#include <stddef.h>

struct prob_cdf {
	double *y;
	size_t len;
};

int prob_cdf_init(struct prob_cdf *cdf, size_t len);

void prob_cdf_free(struct prob_cdf *cdf);

/*
 * Calculate probability distribution (cumulative distribution function) from
 * a set of not normalized, possibly negative samples of length pd->len.
 *
 * The pdf is calculated by calling get_smpl() pd->len times, taking
 * the absolute value, and normalizing to obtain a PDF.
 */
int prob_cdf_from_samples(
	struct prob_cdf *cdf, double (*get_smpl)(void *), void *data);

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
