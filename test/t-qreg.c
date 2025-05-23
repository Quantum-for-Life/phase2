#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "test.h"

#define WD_SEED UINT64_C(0x4809bfb7d258a633)
static struct world WD;

#define WIDTH (64)

#if PHASE2_BACKEND == 0
#define MARGIN (1.0e-14)
#elif PHASE2_BACKEND == 1 /* QuEST */
#define MARGIN (1.0e-11)
#elif PHASE2_BACKEND == 2 /* cuQuantum */
#define MARGIN (1.0e-14)
#endif /* PHASE2_BACKEND */

#define NUM_QUBITS (13)
#define NUM_AMPS (1UL << NUM_QUBITS)
static _Complex double AMPS[NUM_AMPS];

#define SEED UINT64_C(0x34eaaa33)
static struct xoshiro256ss RNG;

static int rand_pauli_op(void)
{
	return (int)(xoshiro256ss_next(&RNG) % 4);
}

static void t_qreg_init(void)
{
	struct qreg reg;
	uint32_t nqb_lo, nqb_hi, wd_size;
	uint64_t namp;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0, "cannot initialize qreg");

	nqb_hi = nqb_lo = 0;
	wd_size = WD.size;
	while (wd_size >>= 1)
		nqb_hi++;
	nqb_lo = NUM_QUBITS - nqb_hi;
	namp = 1UL << nqb_lo;

	TEST_EQ(nqb_lo, reg.qb_lo);
	TEST_EQ(nqb_hi, reg.qb_hi);
	TEST_EQ(namp, reg.namp);

	qreg_free(&reg);
}

static void t_qreg_getsetamp_01(void)
{
	_Complex double z;
	struct qreg reg;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0, "cannot initialize qreg");

	qreg_zero(&reg);
	TEST_ASSERT(reg.namp > 0, "no. of local amps");
	qreg_setamp(&reg, 11, 0.9);
	qreg_setamp(&reg, 111, 0.7 * I);
	qreg_setamp(&reg, 27, -1.9 + 0.8 * I);

	z = 0.0;
	qreg_getamp(&reg, 11, &z);
	TEST_ASSERT(z == 0.9, "amp=11, z=%f+%fi", creal(z), cimag(z));
	z = 0.0;
	qreg_getamp(&reg, 111, &z);
	TEST_ASSERT(z == 0.7 * I, "amp=111, z=%f+%fi", creal(z), cimag(z));
	z = 0.0;
	qreg_getamp(&reg, 27, &z);
	TEST_ASSERT(
		z == -1.9 + 0.8 * I, "amp=27, z=%f+%fi", creal(z), cimag(z));

	qreg_free(&reg);
}

static void t_qreg_getsetamp_02(size_t tag)
{
	_Complex double z;
	struct qreg reg;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0, "cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z = xoshiro256ss_dbl01(&RNG) + xoshiro256ss_dbl01(&RNG) * I;
		AMPS[i] = z;
		qreg_setamp(&reg, i, z);
	}

	for (size_t i = 0; i < NUM_AMPS; i++) {
		qreg_getamp(&reg, i, &z);
		TEST_ASSERT(z == AMPS[i],
			"[%zu] i=%zu, z=%f+%fi, AMPS[i]=%f+%fi", tag, i,
			creal(z), cimag(z), creal(AMPS[i]), cimag(AMPS[i]));
	}

	qreg_free(&reg);
}

static void t_qreg_zero(void)
{
	_Complex double z;
	struct qreg reg;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0, "cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z = xoshiro256ss_dbl01(&RNG) + xoshiro256ss_dbl01(&RNG) * I;
		AMPS[i] = z;
		qreg_setamp(&reg, i, z);
	}

	qreg_zero(&reg);

	for (size_t i = 0; i < NUM_AMPS; i++) {
		qreg_getamp(&reg, i, &z);
		TEST_ASSERT(z == 0.0, "i=%zu, z=%f+%fi", i, creal(z), cimag(z));
	}

	qreg_free(&reg);
}

/*
 * Test if Pauli string equal to identity produces
 * just multiplication by phase.
 */
static void t_qreg_paulirot_00(void)
{
	_Complex double z, z_exp;
	struct qreg reg;
	struct paulis ps_hi, ps_lo;
	double angle;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0, "cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z = xoshiro256ss_dbl01(&RNG) + xoshiro256ss_dbl01(&RNG) * I;
		AMPS[i] = z;
		qreg_setamp(&reg, i, z);
	}

	ps_hi = ps_lo = paulis_new();
	angle = 0.711;
	qreg_paulirot(&reg, ps_hi, &ps_lo, &angle, 1);

	for (size_t i = 0; i < NUM_AMPS; i++) {
		z_exp = AMPS[i] * cexp(I * angle);
		qreg_getamp(&reg, i, &z);
		TEST_ASSERT(cabs(z - z_exp) < MARGIN,
			"i=%zu, z=%f+%fi, z_exp=%f+%fi", i, creal(z), cimag(z),
			creal(z_exp), cimag(z_exp));
	}

	qreg_free(&reg);
}

static void t_qreg_paulirot_0y_explicit(void)
{
	const uint32_t num_qubits = 1;

	_Complex double z, z_exp;
	struct qreg reg;
	struct paulis ps, ps_hi, ps_lo;
	double angle;

	TEST_ASSERT(qreg_init(&reg, num_qubits) == 0, "cannot initialize qreg");
	qreg_zero(&reg);
	qreg_setamp(&reg, 0, 2.0);

	ps = paulis_new();
	paulis_set(&ps, PAULI_Y, 0);
	paulis_split(ps, reg.qb_lo, reg.qb_hi, &ps_lo, &ps_hi);
	angle = 0.11;
	qreg_paulirot(&reg, ps_hi, &ps_lo, &angle, 1);

	z_exp = cos(angle) * 2.0;
	qreg_getamp(&reg, 0, &z);
	TEST_ASSERT(cabs(z - z_exp) < MARGIN, "i=%d, z=%f+%fi, z_exp=%f+%fi", 0,
		creal(z), cimag(z), creal(z_exp), cimag(z_exp));

	z_exp = -sin(angle) * 2.0;
	qreg_getamp(&reg, 1, &z);
	TEST_ASSERT(cabs(z - z_exp) < MARGIN, "i=%d, z=%f+%fi, z_exp=%f+%fi", 1,
		creal(z), cimag(z), creal(z_exp), cimag(z_exp));

	qreg_free(&reg);
}

static void t_qreg_paulirot_0yyy_explicit(void)
{
	const uint32_t num_qubits = 3;

	_Complex double z, z_exp;
	struct qreg reg;
	struct paulis ps, ps_hi, ps_lo;
	double angle;

	TEST_ASSERT(qreg_init(&reg, num_qubits) == 0, "cannot initialize qreg");
	qreg_zero(&reg);
	qreg_setamp(&reg, 0, 2.0);

	ps = paulis_new();
	paulis_set(&ps, PAULI_Y, 0);
	paulis_set(&ps, PAULI_Y, 1);
	paulis_set(&ps, PAULI_Y, 2);
	paulis_split(ps, reg.qb_lo, reg.qb_hi, &ps_lo, &ps_hi);
	angle = 0.11;
	qreg_paulirot(&reg, ps_hi, &ps_lo, &angle, 1);

	z_exp = cos(angle) * 2.0;
	qreg_getamp(&reg, 0, &z);
	TEST_ASSERT(cabs(z - z_exp) < MARGIN, "i=%d, z=%f+%fi, z_exp=%f+%fi", 0,
		creal(z), cimag(z), creal(z_exp), cimag(z_exp));

	z_exp = sin(angle) * 2.0;
	qreg_getamp(&reg, 7, &z);
	TEST_ASSERT(cabs(z - z_exp) < MARGIN, "i=%d, z=%f+%fi, z_exp=%f+%fi", 1,
		creal(z), cimag(z), creal(z_exp), cimag(z_exp));

	qreg_free(&reg);
}

/* Test rotation by one random Pauli string */
static void t_qreg_paulirot_01(size_t tag)
{
	_Complex double z, z_exp;
	struct qreg reg;
	struct paulis ps, ps_hi, ps_lo;
	double angle;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0, "cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z = xoshiro256ss_dbl01(&RNG) + xoshiro256ss_dbl01(&RNG) * I;
		AMPS[i] = z;
		qreg_setamp(&reg, i, z);
	}

	ps = ps_hi = ps_lo = paulis_new();
	for (size_t k = 0; k < NUM_QUBITS; k++)
		paulis_set(&ps, rand_pauli_op(), k);
	paulis_split(ps, reg.qb_lo, reg.qb_hi, &ps_lo, &ps_hi);
	angle = xoshiro256ss_dbl01(&RNG);
	qreg_paulirot(&reg, ps_hi, &ps_lo, &angle, 1);

	for (size_t i = 0; i < NUM_AMPS; i++) {
		_Complex double u = 1.0;

		size_t j = paulis_effect(ps, i, &u);
		z_exp = cos(angle) * AMPS[i] +
			I * conj(u) * sin(angle) * AMPS[j];

		qreg_getamp(&reg, i, &z);
		TEST_ASSERT(cabs(z - z_exp) < MARGIN,
			"[%zu] i=%zu, z=%f+%fi, z_exp=%f+%fi", tag, i, creal(z),
			cimag(z), creal(z_exp), cimag(z_exp));
	}

	qreg_free(&reg);
}

/* Test rotation by two random Pauli strings */
static void t_qreg_paulirot_02(size_t tag)
{
	_Complex double z[2];
	struct qreg reg;
	struct paulis ps[2], ps_hi[2], ps_lo[2];
	double angle[2];

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0, "cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z[0] = xoshiro256ss_dbl01(&RNG) + xoshiro256ss_dbl01(&RNG) * I;
		AMPS[i] = z[0];
		qreg_setamp(&reg, i, z[0]);
	}

	ps[0] = ps_lo[0] = ps_hi[1] = paulis_new();
	ps[1] = ps_lo[1] = ps_hi[1] = paulis_new();
	for (size_t k = 0; k < reg.qb_lo; k++) {
		paulis_set(&ps[0], rand_pauli_op(), k);
		paulis_set(&ps[1], rand_pauli_op(), k);
	}
	for (size_t k = reg.qb_lo; k < reg.qb_lo + reg.qb_hi; k++) {
		int op = rand_pauli_op();
		paulis_set(&ps[0], op, k);
		paulis_set(&ps[1], op, k);
	}
	paulis_split(ps[0], reg.qb_lo, reg.qb_hi, &ps_lo[0], &ps_hi[0]);
	paulis_split(ps[1], reg.qb_lo, reg.qb_hi, &ps_lo[1], &ps_hi[1]);
	TEST_ASSERT(paulis_eq(ps_hi[0], ps_hi[1]),
		"[%zu] hi codes should be equal", tag);

	angle[0] = xoshiro256ss_dbl01(&RNG);
	angle[1] = xoshiro256ss_dbl01(&RNG);
	qreg_paulirot(&reg, ps_hi[0], ps_lo, angle, 2);

	for (size_t k = 0; k < 2; k++) {
		for (size_t i = 0; i < NUM_AMPS; i++) {
			_Complex double u = 1.0;

			size_t j = paulis_effect(ps[k], i, &u);
			if (j < i)
				continue;

			z[0] = cos(angle[k]) * AMPS[i] +
			       I * conj(u) * sin(angle[k]) * AMPS[j];
			z[1] = cos(angle[k]) * AMPS[j] +
			       I * u * sin(angle[k]) * AMPS[i];
			AMPS[i] = z[0];
			AMPS[j] = z[1];
		}
	}

	for (size_t i = 0; i < NUM_AMPS; i++) {
		qreg_getamp(&reg, i, &z[0]);
		TEST_ASSERT(cabs(z[0] - AMPS[i]) < MARGIN,
			"[%zu] i=%zu, z=%f+%fi, z_exp=%f+%fi", tag, i,
			creal(z[0]), cimag(z[0]), creal(AMPS[i]),
			cimag(AMPS[i]));
	}

	qreg_free(&reg);
}

/* Test rotation by n random Pauli strings */
static void t_qreg_paulirot_03(size_t tag, size_t n)
{
	_Complex double z[2];
	struct qreg reg;
	struct paulis *ps, *ps_hi, *ps_lo;
	double *angle;

	ps = malloc(sizeof(struct paulis) * n * 3);
	angle = malloc(sizeof(double) * n);
	if (!ps || !angle) {
		TEST_FAIL("cannot allocate memory");
		return;
	}
	ps_hi = ps + n;
	ps_lo = ps_hi + n;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0, "cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z[0] = xoshiro256ss_dbl01(&RNG) + xoshiro256ss_dbl01(&RNG) * I;
		AMPS[i] = z[0];
		qreg_setamp(&reg, i, z[0]);
	}

	for (size_t l = 0; l < n; l++)
		ps[l] = ps_lo[l] = ps_hi[l] = paulis_new();
	for (size_t k = 0; k < reg.qb_lo; k++)
		for (size_t l = 0; l < n; l++)
			paulis_set(&ps[l], rand_pauli_op(), k);
	for (size_t k = reg.qb_lo; k < reg.qb_lo + reg.qb_hi; k++) {
		int op = rand_pauli_op();
		for (size_t l = 0; l < n; l++)
			paulis_set(&ps[l], op, k);
	}
	for (size_t l = 0; l < n; l++) {
		paulis_split(ps[l], reg.qb_lo, reg.qb_hi, &ps_lo[l], &ps_hi[l]);
		if (l > 0)
			TEST_ASSERT(paulis_eq(ps_hi[0], ps_hi[l]),
				"[%zu] l=%zu hi codes should be equal", tag, l);
	}
	for (size_t l = 0; l < n; l++)
		angle[l] = xoshiro256ss_dbl01(&RNG);
	qreg_paulirot(&reg, ps_hi[0], ps_lo, angle, n);

	for (size_t k = 0; k < n; k++) {
		for (size_t i = 0; i < NUM_AMPS; i++) {
			_Complex double u = 1.0;

			size_t j = paulis_effect(ps[k], i, &u);
			if (j < i)
				continue;

			z[0] = cos(angle[k]) * AMPS[i] +
			       I * conj(u) * sin(angle[k]) * AMPS[j];
			z[1] = cos(angle[k]) * AMPS[j] +
			       I * u * sin(angle[k]) * AMPS[i];
			AMPS[i] = z[0];
			AMPS[j] = z[1];
		}
	}

	for (size_t i = 0; i < NUM_AMPS; i++) {
		qreg_getamp(&reg, i, &z[0]);
		TEST_ASSERT(cabs(z[0] - AMPS[i]) < MARGIN,
			"[%zu] n=%zu, i=%zu, z=%f+%fi, z_exp=%f+%fi", tag, n, i,
			creal(z[0]), cimag(z[0]), creal(AMPS[i]),
			cimag(AMPS[i]));
	}

	qreg_free(&reg);
	free(angle);
	free(ps);
}

static void TEST_MAIN(void)
{
	world_init((void *)0, (void *)0, WD_SEED);
	world_info(&WD);
	log_info("MPI world size: %d", WD.size);

	xoshiro256ss_init(&RNG, SEED);

	t_qreg_init();

	t_qreg_getsetamp_01();
	for (size_t k = 0; k < 99; k++)
		t_qreg_getsetamp_02(k);

	t_qreg_zero();

	t_qreg_paulirot_00();
	if (WD.size == 1) {
		log_info("t_qreg_palirot_0y_explicit()");
		t_qreg_paulirot_0y_explicit();
	}
	if (WD.size <= 4) {
		log_info("t_qreg_palirot_0yyy_explicit()");
		t_qreg_paulirot_0yyy_explicit();
	}
	t_qreg_paulirot_01(13);

	for (size_t k = 0; k < 99; k++)
		t_qreg_paulirot_01(k);
	for (size_t k = 0; k < 99; k++)
		t_qreg_paulirot_02(k);
	for (size_t k = 0; k < 9; k++)
		for (size_t n = 1; n <= 99; n++)
			t_qreg_paulirot_03(k, n);
	t_qreg_paulirot_03(1234, 999);

	world_free();
}
