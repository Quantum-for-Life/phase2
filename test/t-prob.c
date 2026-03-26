#include "c23_compat.h"
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "phase2/prob.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0xb3f7a28e510cd946)

struct iter_data {
	size_t i;
	double *vals;
};

static double iter_cb(void *data)
{
	struct iter_data *d = data;
	return d->vals[d->i++];
}

/*
 * 4-element uniform distribution [1, 1, 1, 1].
 * CDF: [0.25, 0.5, 0.75, 1.0].
 */
static void t_prob_cdf_uniform(void)
{
	struct prob_cdf cdf;
	double vals[] = { 1.0, 1.0, 1.0, 1.0 };
	struct iter_data d = { .i = 0, .vals = vals };

	TEST_EQ(prob_cdf_init(&cdf, 4), 0);
	TEST_EQ(prob_cdf_from_iter(&cdf, iter_cb, &d), 0);

	/* CDF must be monotonically non-decreasing. */
	for (size_t k = 1; k < cdf.len; k++)
		TEST_ASSERT(cdf.y[k] >= cdf.y[k - 1],
			"CDF not monotone at k=%zu", k);

	/* Last element must be 1.0 (within floating-point tolerance). */
	TEST_ASSERT(fabs(cdf.y[3] - 1.0) < DBL_EPSILON,
		"CDF last element %f != 1.0", cdf.y[3]);

	/* Inverse: y=0 maps to index 0. */
	TEST_EQ(prob_cdf_inverse(&cdf, 0.0), (size_t)0);

	/* Inverse: y=1.0 maps to last index. */
	TEST_EQ(prob_cdf_inverse(&cdf, 1.0), (size_t)3);

	/* Inverse: y just above 0.25 maps past index 0. */
	TEST_EQ(prob_cdf_inverse(&cdf, 0.26), (size_t)1);

	/* Inverse: y just above 0.5 maps past index 1. */
	TEST_EQ(prob_cdf_inverse(&cdf, 0.51), (size_t)2);

	prob_cdf_free(&cdf);
}

/*
 * Single-element CDF.  Inverse must always return 0.
 */
static void t_prob_cdf_single(void)
{
	struct prob_cdf cdf;
	double vals[] = { 5.0 };
	struct iter_data d = { .i = 0, .vals = vals };

	TEST_EQ(prob_cdf_init(&cdf, 1), 0);
	TEST_EQ(prob_cdf_from_iter(&cdf, iter_cb, &d), 0);

	TEST_ASSERT(fabs(cdf.y[0] - 1.0) < DBL_EPSILON,
		"single-element CDF should be 1.0");

	TEST_EQ(prob_cdf_inverse(&cdf, 0.0), (size_t)0);
	TEST_EQ(prob_cdf_inverse(&cdf, 0.5), (size_t)0);
	TEST_EQ(prob_cdf_inverse(&cdf, 1.0), (size_t)0);

	prob_cdf_free(&cdf);
}

/*
 * 3 elements with weights [1, 2, 1].  Sum = 4.
 * PDF: [0.25, 0.5, 0.25].
 * CDF: [0.25, 0.75, 1.0].
 */
static void t_prob_cdf_known(void)
{
	struct prob_cdf cdf;
	double vals[] = { 1.0, 2.0, 1.0 };
	struct iter_data d = { .i = 0, .vals = vals };

	TEST_EQ(prob_cdf_init(&cdf, 3), 0);
	TEST_EQ(prob_cdf_from_iter(&cdf, iter_cb, &d), 0);

	TEST_ASSERT(fabs(cdf.y[0] - 0.25) < 1e-12,
		"CDF[0] = %f, expected 0.25", cdf.y[0]);
	TEST_ASSERT(fabs(cdf.y[1] - 0.75) < 1e-12,
		"CDF[1] = %f, expected 0.75", cdf.y[1]);
	TEST_ASSERT(fabs(cdf.y[2] - 1.0) < 1e-12,
		"CDF[2] = %f, expected 1.0", cdf.y[2]);

	/* y=0.0: below all CDF values, returns 0. */
	TEST_EQ(prob_cdf_inverse(&cdf, 0.0), (size_t)0);

	/*
	 * y=0.25: equal to CDF[0], the code advances past it,
	 * returns 1.
	 */
	TEST_EQ(prob_cdf_inverse(&cdf, 0.25), (size_t)1);

	/* y=0.5: between CDF[0] and CDF[1], returns 1. */
	TEST_EQ(prob_cdf_inverse(&cdf, 0.5), (size_t)1);

	/* y=0.75: equal to CDF[1], advances past it, returns 2. */
	TEST_EQ(prob_cdf_inverse(&cdf, 0.75), (size_t)2);

	/* y=1.0: last index. */
	TEST_EQ(prob_cdf_inverse(&cdf, 1.0), (size_t)2);

	prob_cdf_free(&cdf);
}

/*
 * Negative weights [-1, 2, -1] produce the same CDF as [1, 2, 1]
 * because absolute values are taken.
 */
static void t_prob_cdf_negative(void)
{
	struct prob_cdf cdf;
	double vals[] = { -1.0, 2.0, -1.0 };
	struct iter_data d = { .i = 0, .vals = vals };

	TEST_EQ(prob_cdf_init(&cdf, 3), 0);
	TEST_EQ(prob_cdf_from_iter(&cdf, iter_cb, &d), 0);

	TEST_ASSERT(fabs(cdf.y[0] - 0.25) < 1e-12,
		"CDF[0] = %f, expected 0.25", cdf.y[0]);
	TEST_ASSERT(fabs(cdf.y[1] - 0.75) < 1e-12,
		"CDF[1] = %f, expected 0.75", cdf.y[1]);
	TEST_ASSERT(fabs(cdf.y[2] - 1.0) < 1e-12,
		"CDF[2] = %f, expected 1.0", cdf.y[2]);

	prob_cdf_free(&cdf);
}

/*
 * All-zero weights: prob_cdf_from_iter must return -1.
 */
static void t_prob_cdf_zeros(void)
{
	struct prob_cdf cdf;
	double vals[] = { 0.0, 0.0, 0.0 };
	struct iter_data d = { .i = 0, .vals = vals };

	TEST_EQ(prob_cdf_init(&cdf, 3), 0);
	TEST_EQ(prob_cdf_from_iter(&cdf, iter_cb, &d), -1);

	prob_cdf_free(&cdf);
}

/*
 * Boundary tests for prob_cdf_inverse at y=0 and y=1.
 */
static void t_prob_cdf_inverse_boundaries(void)
{
	struct prob_cdf cdf;
	double vals[] = { 1.0, 3.0, 1.0 };
	struct iter_data d = { .i = 0, .vals = vals };

	/* CDF: [0.2, 0.8, 1.0] */
	TEST_EQ(prob_cdf_init(&cdf, 3), 0);
	TEST_EQ(prob_cdf_from_iter(&cdf, iter_cb, &d), 0);

	/* y=0: strictly below CDF[0]=0.2, returns 0. */
	TEST_EQ(prob_cdf_inverse(&cdf, 0.0), (size_t)0);

	/* y=1.0: at or above all CDF values, returns last index. */
	TEST_EQ(prob_cdf_inverse(&cdf, 1.0), (size_t)2);

	/* y slightly above 0 still returns 0. */
	TEST_EQ(prob_cdf_inverse(&cdf, 1e-15), (size_t)0);

	/* y just below 1.0: CDF[1]=0.8 <= 0.999, advance to 2. */
	TEST_EQ(prob_cdf_inverse(&cdf, 0.999), (size_t)2);

	prob_cdf_free(&cdf);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	t_prob_cdf_uniform();
	t_prob_cdf_single();
	t_prob_cdf_known();
	t_prob_cdf_negative();
	t_prob_cdf_zeros();
	t_prob_cdf_inverse_boundaries();

	world_free();
}
