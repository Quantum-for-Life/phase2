#include "c23_compat.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "det_small.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "test.h"

#define WD_SEED UINT64_C(0x7a3f88fc97e2c811)
#define RNG_SEED UINT64_C(0xc1da4eb88aa72cd1)

#define MAX_N (17u)

static struct xoshiro256ss RNG;

static double rand_uniform(void)
{
	return 2.0 * xoshiro256ss_dbl01(&RNG) - 1.0;
}

static double det_naive(const double *A, uint32_t n)
{
	if (n == 0)
		return 1.0;
	if (n == 1)
		return A[0];
	if (n == 2)
		return A[0] * A[3] - A[1] * A[2];

	double M[MAX_N * MAX_N];
	for (uint32_t i = 0; i < n * n; i++)
		M[i] = A[i];

	double det = 1.0;
	for (uint32_t k = 0; k < n; k++) {
		uint32_t piv = k;
		double mx = fabs(M[k * n + k]);
		for (uint32_t i = k + 1; i < n; i++)
			if (fabs(M[i * n + k]) > mx) {
				mx = fabs(M[i * n + k]);
				piv = i;
			}
		if (mx == 0.0)
			return 0.0;
		if (piv != k) {
			for (uint32_t j = k; j < n; j++) {
				double t = M[k * n + j];
				M[k * n + j] = M[piv * n + j];
				M[piv * n + j] = t;
			}
			det = -det;
		}
		double p = M[k * n + k];
		det *= p;
		for (uint32_t i = k + 1; i < n; i++) {
			double f = M[i * n + k] / p;
			for (uint32_t j = k + 1; j < n; j++)
				M[i * n + j] -= f * M[k * n + j];
		}
	}
	return det;
}

static void t_edge_cases(void)
{
	TEST_ASSERT(det_small(NULL, 0) == 1.0, "n=0 det");

	double A1[1] = { 3.5 };
	TEST_ASSERT(det_small(A1, 1) == 3.5, "n=1 det");

	double Id[3 * 3] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
	TEST_ASSERT(det_small(Id, 3) == 1.0, "n=3 identity");

	double Z[4 * 4] = { 0 };
	TEST_ASSERT(det_small(Z, 4) == 0.0, "n=4 zero");

	double P[3 * 3] = { 0, 1, 0, 0, 0, 1, 1, 0, 0 };
	const double dP = det_small(P, 3);
	TEST_ASSERT(dP == 1.0 || dP == -1.0, "permutation det = +/-1, got %f",
		dP);

	double U[4 * 4] = {
		2.0, 1.0, 3.0, 4.0,
		0.0, -1.0, 2.0, 1.0,
		0.0, 0.0, 5.0, 2.0,
		0.0, 0.0, 0.0, -3.0
	};
	const double dU = det_small(U, 4);
	const double exp = 2.0 * -1.0 * 5.0 * -3.0;
	TEST_NEAR(dU, exp, 1e-12);
}

static void t_random(uint32_t n, size_t trials)
{
	const size_t sz = (size_t)n * n;
	double *A = malloc(sizeof(double) * sz);
	for (size_t t = 0; t < trials; t++) {
		for (size_t i = 0; i < sz; i++)
			A[i] = rand_uniform();
		const double d_ref = det_naive(A, n);
		const double d_got = det_small(A, n);
		const double scale = fabs(d_ref) > 1.0 ? fabs(d_ref) : 1.0;
		TEST_ASSERT(fabs(d_got - d_ref) < 1e-10 * scale,
			"n=%u trial=%zu got=%g ref=%g", n, t, d_got, d_ref);
	}
	free(A);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	xoshiro256ss_init(&RNG, RNG_SEED);

	t_edge_cases();
	for (uint32_t n = 1; n <= MAX_N; n++)
		t_random(n, 32);

	world_free();
}
