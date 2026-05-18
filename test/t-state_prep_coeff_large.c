/*
 * Slow opt-in test: N=14 closed-shell expansion stress.  This
 * is gated behind `make test-slow` and not part of the
 * default check target.
 */

#define _POSIX_C_SOURCE 200809L
#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "phase2/qreg.h"
#include "phase2/state_prep_coeff.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "test.h"

#define WD_SEED UINT64_C(0xaa881e6c992b7411)
#define RNG_SEED UINT64_C(0xc77f1109aa8de822)

#define N_SITES 14u
#define N_OCC 7u

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	struct xoshiro256ss rng;
	xoshiro256ss_init(&rng, RNG_SEED);

	double Ca[N_SITES * N_OCC];
	for (uint32_t i = 0; i < N_SITES * N_OCC; i++)
		Ca[i] = 2.0 * xoshiro256ss_dbl01(&rng) - 1.0;
	/* Crude column-orthonormalisation by Gram-Schmidt. */
	for (uint32_t c = 0; c < N_OCC; c++) {
		for (uint32_t c2 = 0; c2 < c; c2++) {
			double dot = 0.0;
			for (uint32_t r = 0; r < N_SITES; r++)
				dot += Ca[r * N_OCC + c]
				       * Ca[r * N_OCC + c2];
			for (uint32_t r = 0; r < N_SITES; r++)
				Ca[r * N_OCC + c] -=
					dot * Ca[r * N_OCC + c2];
		}
		double norm = 0.0;
		for (uint32_t r = 0; r < N_SITES; r++)
			norm += Ca[r * N_OCC + c] * Ca[r * N_OCC + c];
		norm = sqrt(norm);
		for (uint32_t r = 0; r < N_SITES; r++)
			Ca[r * N_OCC + c] /= norm;
	}

	struct qreg reg;
	TEST_EQ(qreg_init(&reg, 2 * N_SITES), 0);
	qreg_zero(&reg);

	struct state_prep_coeff_scratch sc;
	TEST_EQ(state_prep_coeff_scratch_init(&sc, N_SITES, N_OCC, N_OCC), 0);

	struct timespec t1, t2;
	clock_gettime(CLOCK_MONOTONIC, &t1);
	TEST_EQ(state_prep_coeff_expand(&reg, &sc, Ca, NULL, 1.0, 0, 0), 0);
	clock_gettime(CLOCK_MONOTONIC, &t2);
	const double dt =
		(t2.tv_sec - t1.tv_sec)
		+ (t2.tv_nsec - t1.tv_nsec) * 1e-9;
	fprintf(stderr, "N=14 expand: %.3f s\n", dt);
	TEST_ASSERT(dt < 60.0, "expansion exceeded 60s budget: %.3f", dt);

	/*
	 * Sanity check at scale: the trial state must be
	 * normalised to within a loose 1e-10 tolerance.  The
	 * tolerance absorbs FP accumulation across ~2.6 M terms;
	 * a stricter bound would be optimistic.  We use
	 * state_prep_coeff_inner with the same C to compute
	 * <psi|psi> directly, which is independent of the
	 * scatter-and-gather path qreg_getamp would take.
	 */
	_Complex double inner = 0.0;
	TEST_EQ(state_prep_coeff_inner(&reg, &sc, Ca, NULL, 1.0, 0, &inner), 0);
	state_prep_coeff_scratch_free(&sc);
	const double nsq = creal(inner);
	fprintf(stderr, "N=14 <psi|psi>: %.6f (imag=%.2e)\n", nsq,
		cimag(inner));
	TEST_NEAR(nsq, 1.0, 1e-10);
	TEST_NEAR(cimag(inner), 0.0, 1e-10);

	qreg_free(&reg);
	world_free();
}
