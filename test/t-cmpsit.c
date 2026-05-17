/*
 * t-cmpsit -- correctness checks for the composite
 * (partially-randomised) 2nd-order Trotter
 * (circ/cmpsit.c).
 *
 * Single scenario today: seed reproducibility.  Two
 * cmpsit_simul runs with identical cmpsit_data.seed must
 * produce bit-identical per-sample overlaps.  Pins the
 * RNG seed contract -- the same seed must drive the same
 * sample sequence across all `samples` draws.
 *
 * Mean-convergence checks against a trotter / trott2
 * reference are the natural next step once a clean
 * oracle for the partial-randomisation case lands.  (A
 * fully-deterministic cmpsit -- length = hm.len, depth =
 * 0 -- would be the obvious degenerate case, but it
 * trips prob_cdf_init's len > 0 guard.)
 *
 * Fixtures programmatic; no HDF5.
 */

#define LOG_SUBSYS "test"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "circ/cmpsit.h"
#include "log.h"
#include "phase2.h"
#include "xoshiro256ss.h"

#include "test.h"

#define WD_SEED   UINT64_C(0x73a82de0c1f59048)
#define RNG_SEED  UINT64_C(0x4f60c9837ea51d62)

#define NQB       5
#define NUM_TERMS 6
#define NUM_DETS  4

static struct xoshiro256ss RNG;


/* -- programmatic fixture builders ----------------------------------------*/

static struct circ_hamil build_hamil(size_t nterms)
{
	struct circ_hamil hm;
	TEST_EQ(circ_hamil_init(&hm, NQB, nterms), 0);
	for (size_t i = 0; i < nterms; i++) {
		hm.terms[i].op = paulis_new();
		for (uint32_t k = 0; k < NQB; k++)
			paulis_set(&hm.terms[i].op,
				(enum pauli_op)(xoshiro256ss_next(&RNG) % 4),
				k);
		hm.terms[i].cf = xoshiro256ss_dbl01(&RNG) - 0.5;
	}
	return hm;
}

static struct circ_muldet build_muldet(size_t ndets)
{
	struct circ_muldet md;
	TEST_EQ(circ_muldet_init(&md, ndets), 0);
	const uint64_t mask = (NQB >= 64)
		? ~UINT64_C(0)
		: ((UINT64_C(1) << NQB) - 1);
	double norm = 0.0;
	for (size_t i = 0; i < ndets; i++) {
		md.dets[i].idx = xoshiro256ss_next(&RNG) & mask;
		const double re = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;
		const double im = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;
		md.dets[i].cf = re + im * I;
		norm += creal(md.dets[i].cf) * creal(md.dets[i].cf) +
			cimag(md.dets[i].cf) * cimag(md.dets[i].cf);
	}
	norm = sqrt(norm);
	for (size_t i = 0; i < ndets; i++)
		md.dets[i].cf /= norm;
	return md;
}

static struct circ_hamil clone_hamil(const struct circ_hamil *src)
{
	struct circ_hamil copy;
	TEST_EQ(circ_hamil_init(&copy, src->qb, src->len), 0);
	for (size_t i = 0; i < src->len; i++)
		copy.terms[i] = src->terms[i];
	return copy;
}

static struct circ_muldet clone_muldet(const struct circ_muldet *src)
{
	struct circ_muldet copy;
	TEST_EQ(circ_muldet_init(&copy, src->len), 0);
	for (size_t i = 0; i < src->len; i++)
		copy.dets[i] = src->dets[i];
	return copy;
}


/* -- tests -----------------------------------------------------------------*/

static void t_seed_reproducible(void)
{
	struct circ_hamil hm_master = build_hamil(NUM_TERMS);
	struct circ_muldet md_master = build_muldet(NUM_DETS);

	const struct cmpsit_data dt = {
		.samples = 3,
		.steps = 2,
		.length = 2,
		.depth = 4,
		.angle_det = 0.1,
		.angle_rand = 0.1,
		.seed = UINT64_C(0xc6a51e3b80f72d94),
	};

	struct cmpsit a, b;
	struct circ_hamil hm_a = clone_hamil(&hm_master);
	struct circ_muldet md_a = clone_muldet(&md_master);
	TEST_EQ(cmpsit_init(&a, &dt, hm_a, STPREP_MULTIDET, &md_a, NULL),
		0);
	struct circ_hamil hm_b = clone_hamil(&hm_master);
	struct circ_muldet md_b = clone_muldet(&md_master);
	TEST_EQ(cmpsit_init(&b, &dt, hm_b, STPREP_MULTIDET, &md_b, NULL),
		0);

	TEST_EQ(cmpsit_simul(&a), 0);
	TEST_EQ(cmpsit_simul(&b), 0);

	for (size_t i = 0; i < a.ct.vals.len; i++) {
		const _Complex double za = a.ct.vals.z[i];
		const _Complex double zb = b.ct.vals.z[i];
		TEST_ASSERT(za == zb,
			"sample %zu diverges: a=(%g+%gi) b=(%g+%gi)",
			i, creal(za), cimag(za), creal(zb), cimag(zb));
	}

	cmpsit_free(&a);
	cmpsit_free(&b);

	circ_hamil_free(&hm_master);
	circ_muldet_free(&md_master);
}


int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	xoshiro256ss_init(&RNG, RNG_SEED);

	t_seed_reproducible();

	world_free();
	return 0;
}
