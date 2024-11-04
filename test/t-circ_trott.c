#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "circ/trott.h"
#include "phase2.h"
#include "xoshiro256ss.h"

#include "test.h"

#define WD_SEED UINT64_C(0x682011f6dd97fc67)
static struct world WD;

#if PHASE2_BACKEND == 0
#define MARGIN (1.0e-14)
#elif PHASE2_BACKEND == 1 /* QuEST */
#define MARGIN (1.0e-6)
#elif PHASE2_BACKEND == 2 /* cuQuantum */
#define MARGIN (1.0e-14)
#endif /* PHASE2_BACKEND */

#define WIDTH (64)
#define NUM_QUBITS (12)
#define NUM_AMPS (1UL << NUM_QUBITS)
static _Complex double AMPS_INIT[NUM_AMPS];
static _Complex double AMPS[NUM_AMPS];

#define TROTT_STEPS_MAX (99UL)
static size_t TROTT_STEPS;
static _Complex double TROTT_VALS[TROTT_STEPS_MAX];

#define HAMIL_TERMS_MAX (99UL)
static size_t HAMIL_TERMS;
static double HAMIL_COEFFS[HAMIL_TERMS_MAX];
static struct paulis HAMIL_PAULIS[HAMIL_TERMS_MAX];
static double HAMIL_DELTA = 1.0;

#define MULTIDET_DETS_MAX (99UL)
static size_t MULTIDET_DETS;
static _Complex double MULTIDET_COEFFS[MULTIDET_DETS_MAX];
static size_t MULTIDET_IDX[MULTIDET_DETS_MAX];

#define SEED UINT64_C(0xd3b9268b8737ddc0)
static struct xoshiro256ss RNG;

extern int trott_simulate(struct circ *c);

static double rand_double(void)
{
	uint64_t x = xoshiro256ss_next(&RNG);

	return (x >> 11) * 0x1.0p-53;
}

static _Complex double rand_complex(void)
{
	_Complex double z = CMPLX(rand_double(), rand_double());

	return z / cabs(z);
}

static int rand_pauli(void)
{
	int x = (int)(xoshiro256ss_next(&RNG) % 4);

	return x;
}

static struct paulis rand_paulis(uint32_t num_qubits)
{
	struct paulis ps = paulis_new();
	for (uint32_t k = 0; k < num_qubits; k++)
		paulis_set(&ps, rand_pauli(), k);

	return ps;
}

/* Mockup trotter circuit. Explicit and local. */
static void trotter_mockup(void)
{
	for (size_t i = 0; i < NUM_AMPS; i++)
		AMPS[i] = AMPS_INIT[i];

	for (size_t s = 0; s < TROTT_STEPS; s++) {
		for (size_t k = 0; k < HAMIL_TERMS; k++) {
			for (size_t i = 0; i < NUM_AMPS; i++) {
				_Complex double z = 1.0;
				size_t j =
					paulis_effect(HAMIL_PAULIS[k], i, &z);
				if (j < i)
					continue;

				double ph = HAMIL_COEFFS[k] * HAMIL_DELTA;
				_Complex double x, y;
				x = cos(ph) * AMPS[i] +
				    I * conj(z) * sin(ph) * AMPS[j];
				y = cos(ph) * AMPS[j] +
				    I * z * sin(ph) * AMPS[i];

				AMPS[i] = x;
				AMPS[j] = y;
			}
		}

		TROTT_VALS[s] = 0.0;
		for (size_t i = 0; i < NUM_AMPS; i++)
			TROTT_VALS[s] += conj(AMPS_INIT[i]) * AMPS[i];
	}
}

static void t_trotter_mockup_sanity_check(void)
{
	AMPS_INIT[0] = 1.0;
	for (size_t i = 1; i < NUM_AMPS; i++)
		AMPS_INIT[i] = 0.0;

	TROTT_STEPS = 1;

	HAMIL_DELTA = 1.0;
	HAMIL_TERMS = 1;
	HAMIL_PAULIS[0] = paulis_new();
	HAMIL_COEFFS[0] = 1.0;

	trotter_mockup();

	TEST_ASSERT(TROTT_VALS[0] == cexp(CMPLX(0.0, 1.0)),
		"trotter circuit sanity check");
}

static void t_circ_trott(size_t tag, size_t ts, size_t md, size_t ht)
{
	TROTT_STEPS = ts;
	MULTIDET_DETS = md;
	HAMIL_TERMS = ht;
	HAMIL_DELTA = rand_double() * 0.5 + 0.5;

	/* Generate random multidet / hamiltonian */
	double norm = 0.0;
	for (size_t m = 0; m < MULTIDET_DETS; m++) {
		MULTIDET_COEFFS[m] = rand_complex();
		double x = cabs(MULTIDET_COEFFS[m]);
		norm += x * x;
		MULTIDET_IDX[m] = m;
	}
	norm = sqrt(norm);
	for (size_t m = 0; m < MULTIDET_DETS; m++)
		MULTIDET_COEFFS[m] /= norm;
	for (size_t k = 0; k < HAMIL_TERMS; k++) {
		HAMIL_COEFFS[k] = rand_double();
		HAMIL_PAULIS[k] = rand_paulis(NUM_QUBITS);
	}

	/* Initialize state for mockup circuit. */
	for (size_t i = 0; i < NUM_AMPS; i++)
		AMPS_INIT[i] = 0.0;
	for (size_t m = 0; m < MULTIDET_DETS; m++)
		AMPS_INIT[MULTIDET_IDX[m]] = MULTIDET_COEFFS[m];

	trotter_mockup();

	/* Initialize circ_trott manually, because we don't have a
	 * handle to an open data file. */
	struct circ_trott tt;
	tt.circ.simulate = trott_simulate;
	qreg_init(&tt.circ.reg, NUM_QUBITS);
	circ_cache_init(&tt.circ.cache, tt.circ.reg.qb_lo, tt.circ.reg.qb_hi);
	tt.delta = HAMIL_DELTA;

	struct circ_hamil *h = &tt.circ.hamil;
	h->nqb = NUM_QUBITS;
	h->nterms = HAMIL_TERMS;
	h->terms = malloc(sizeof *h->terms * HAMIL_TERMS);
	if (!h->terms) {
		TEST_FAIL("malloc h->terms");
		goto malloc_hterms;
	}
	for (size_t k = 0; k < HAMIL_TERMS; k++) {
		h->terms[k].cf = HAMIL_COEFFS[k];
		h->terms[k].op = HAMIL_PAULIS[k];
	}

	struct circ_muldet *m = &tt.circ.muldet;
	m->ndets = MULTIDET_DETS;
	m->dets = malloc(sizeof *m->dets * MULTIDET_DETS);
	if (!m->dets) {
		TEST_FAIL("malloc m->dets");
		goto malloc_mddets;
	}
	for (size_t k = 0; k < MULTIDET_DETS; k++) {
		m->dets[k].idx = MULTIDET_IDX[k];
		m->dets[k].cf = MULTIDET_COEFFS[k];
	}

	tt.res.nsteps = TROTT_STEPS;
	tt.res.steps = malloc(sizeof *tt.res.steps * TROTT_STEPS);
	if (!tt.res.steps) {
		TEST_FAIL("malloc res.steps");
		goto malloc_ressteps;
	}

	circ_simulate(&tt.circ);

	/* Compare results. */
	for (size_t s = 0; s < TROTT_STEPS; s++) {
		_Complex double eval = tt.res.steps[s];
		TEST_ASSERT(cabs(tt.res.steps[s] - TROTT_VALS[s]) < MARGIN,
			"[%zu] ts=%zu, md=%zu, ht=%zu, step=%zu "
			"eval=%f+%fi, TROTT_VALS=%f+%fi",
			tag, ts, md, ht, s, creal(eval), cimag(eval),
			creal(TROTT_VALS[s]), cimag(TROTT_VALS[s]));
	}

	// free(res.steps);  -- freed by circ_destroy()
malloc_ressteps:
	free(m->dets);
malloc_mddets:
	free(h->terms);
malloc_hterms:;
}

void TEST_MAIN(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&WD);
	log_info("MPI world size: %d", WD.size);

	xoshiro256ss_init(&RNG, SEED);

	t_trotter_mockup_sanity_check();

	for (size_t n = 0; n < 9; n++)
		for (size_t ts = 1; ts < TROTT_STEPS_MAX; ts++) {
			t_circ_trott(n, ts, 1, 1);
			t_circ_trott(n, ts, 5, 1);
			t_circ_trott(n, ts, 1, 5);
		}

	for (size_t n = 0; n < 9; n++)
		for (size_t ht = 1; ht < HAMIL_TERMS_MAX; ht++) {
			t_circ_trott(n, 1, 1, ht);
			t_circ_trott(n, 7, 1, ht);
			t_circ_trott(n, 1, 7, ht);
		}

	for (size_t n = 0; n < 9; n++)
		for (size_t md = 1; md < MULTIDET_DETS_MAX; md++) {
			t_circ_trott(n, 1, md, 1);
			t_circ_trott(n, 3, md, 1);
			t_circ_trott(n, 1, md, 3);
		}

	t_circ_trott(9999, TROTT_STEPS_MAX - 1, MULTIDET_DETS_MAX - 1,
		HAMIL_TERMS_MAX - 1);

	world_destroy();
}
