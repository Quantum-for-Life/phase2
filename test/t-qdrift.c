/*
 * t-qdrift -- correctness checks for the qDRIFT
 * randomised product formula (circ/qdrift.c).
 *
 * Two scenarios:
 *
 *   t_seed_reproducible    Two qdrift_simul runs with
 *                          identical seeds must produce
 *                          bit-identical per-sample
 *                          overlaps.  Pins the RNG seed
 *                          contract.
 *
 *   t_single_term_eq_trott Single-term Hamiltonian: every
 *                          qdrift draw selects the same
 *                          term, so the sampler degenerates
 *                          to a deterministic Trotter sweep.
 *                          The per-sample overlap must
 *                          match a first-order Trotter run
 *                          at matching total angle to
 *                          floating-point precision.
 *
 * Fixtures are built in-process; no HDF5 dependency.
 */

#define LOG_SUBSYS "test"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "circ/qdrift.h"
#include "circ/trott.h"
#include "log.h"
#include "phase2.h"
#include "xoshiro256ss.h"

#include "test.h"

#define WD_SEED   UINT64_C(0xe6b81f0c4a39271d)
#define RNG_SEED  UINT64_C(0x927b4ce30a8df615)

#define NQB         5
#define NUM_TERMS   6
#define NUM_DETS    4

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

/* Allocate a deep copy so the second qdrift_init has independent
 * buffers.  circ takes ownership; both halves end up freed. */
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
	/* Build a master fixture; clone for each qdrift instance
	 * since circ_init adopts ownership. */
	struct circ_hamil hm_master = build_hamil(NUM_TERMS);
	struct circ_muldet md_master = build_muldet(NUM_DETS);

	const struct qdrift_data dt = {
		.samples = 4,
		.depth = 8,
		.step_size = 0.1,
		.seed = UINT64_C(0xa5f72e0d8b41069c),
	};

	struct qdrift a, b;
	struct circ_hamil hm_a = clone_hamil(&hm_master);
	struct circ_muldet md_a = clone_muldet(&md_master);
	TEST_EQ(qdrift_init(&a, &dt, hm_a, STPREP_MULTIDET, &md_a, NULL),
		0);
	struct circ_hamil hm_b = clone_hamil(&hm_master);
	struct circ_muldet md_b = clone_muldet(&md_master);
	TEST_EQ(qdrift_init(&b, &dt, hm_b, STPREP_MULTIDET, &md_b, NULL),
		0);

	TEST_EQ(qdrift_simul(&a), 0);
	TEST_EQ(qdrift_simul(&b), 0);

	for (size_t i = 0; i < a.ct.vals.len; i++) {
		const _Complex double za = a.ct.vals.z[i];
		const _Complex double zb = b.ct.vals.z[i];
		TEST_ASSERT(za == zb,
			"sample %zu diverges: a=(%g+%gi) b=(%g+%gi)",
			i, creal(za), cimag(za), creal(zb), cimag(zb));
	}

	qdrift_free(&a);
	qdrift_free(&b);

	circ_hamil_free(&hm_master);
	circ_muldet_free(&md_master);
}

static void t_single_term_eq_trott(void)
{
	/* Single-term Hamiltonian H = +1 * P.  qdrift always
	 * samples this term, carrying signof(c0) = +1 onto the
	 * random term.  A depth-deep qdrift sample then applies
	 * exp(-i * asin(step_size) * P) `depth` times.
	 *
	 * trott with steps=1, delta = depth * asin(step_size)
	 * applies exp(-i * delta * c0 * P) = exp(-i * total
	 * * P) -- identical to the qdrift sweep since c0 = 1
	 * and Pauli powers compose without further error.
	 *
	 * Compare the per-sample overlap. */
	const size_t depth = 4;
	const double step_size = 0.2;
	const double delta = depth * asin(step_size);

	/* Hand-built single-term Hamiltonian with c0 = 1 so the
	 * sign carried by qdrift matches the magnitude trott
	 * applies. */
	struct circ_hamil hm_q;
	TEST_EQ(circ_hamil_init(&hm_q, NQB, 1), 0);
	hm_q.terms[0].op = paulis_new();
	for (uint32_t k = 0; k < NQB; k++)
		paulis_set(&hm_q.terms[0].op,
			(enum pauli_op)(xoshiro256ss_next(&RNG) % 4), k);
	hm_q.terms[0].cf = 1.0;

	struct circ_muldet md_q = build_muldet(NUM_DETS);
	struct circ_hamil hm_t = clone_hamil(&hm_q);
	struct circ_muldet md_t = clone_muldet(&md_q);

	const struct qdrift_data qdt = {
		.samples = 1,
		.depth = depth,
		.step_size = step_size,
		.seed = UINT64_C(0xc9f30b51d7204e6a),
	};
	struct qdrift q;
	TEST_EQ(qdrift_init(&q, &qdt, hm_q, STPREP_MULTIDET, &md_q, NULL),
		0);
	TEST_EQ(qdrift_simul(&q), 0);

	const struct trott_data tdt = {
		.steps = 1,
		.delta = delta,
	};
	struct trott t;
	TEST_EQ(trott_init(&t, &tdt, hm_t, STPREP_MULTIDET, &md_t, NULL),
		0);
	TEST_EQ(trott_simul(&t), 0);

	const _Complex double zq = q.ct.vals.z[0];
	const _Complex double zt = t.ct.vals.z[0];
	TEST_ASSERT(cabs(zq - zt) < 1e-12,
		"qdrift /= trott on single-term: "
		"q=(%g+%gi) t=(%g+%gi) |q-t|=%g",
		creal(zq), cimag(zq), creal(zt), cimag(zt),
		cabs(zq - zt));

	qdrift_free(&q);
	trott_free(&t);
}


int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	xoshiro256ss_init(&RNG, RNG_SEED);

	t_seed_reproducible();
	t_single_term_eq_trott();

	world_free();
	return 0;
}
