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

static struct world WD;

#define WIDTH (64)
#define MARGIN (1.0e-10)

#define NUM_QUBITS (7)
#define NUM_AMPS (1UL << NUM_QUBITS)
static _Complex double AMPS[NUM_AMPS];
#define NUM_PAULIS (99)
static struct paulis PS[NUM_PAULIS];
static double angles[NUM_PAULIS];


#define SEED UINT64_C(0x34eaaa33)
static struct xoshiro256ss RNG;

double rand_dbl(void)
{
	uint64_t x = xoshiro256ss_next(&RNG);

	return (x >> 11) * 0x1.0p-53;
}

enum pauli_op rand_pauli_op(void)
{
	enum pauli_op x = (int)(xoshiro256ss_next(&RNG) % 4);
	
	return x;
}

void test_qreg_init(void)
{
	struct qreg reg;
	uint32_t qb_lo, qb_hi, wd_size;
	uint64_t num_amps;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0,
		"cannot initialize qreg");

	qb_hi = qb_lo = 0;
	wd_size = WD.size;
	while (wd_size >>= 1)
		qb_hi++;
	qb_lo = NUM_QUBITS - qb_hi;
	num_amps = 1UL << qb_lo;

	TEST_EQ(qb_lo, reg.qb_lo);
	TEST_EQ(qb_hi, reg.qb_hi);
	TEST_EQ(num_amps, reg.num_amps);

	qreg_destroy(&reg);
}

void test_qreg_getsetamp_01(void)
{
	_Complex double z;
	struct qreg reg;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0,
		"cannot initialize qreg");

	qreg_zero(&reg);
	TEST_ASSERT(reg.num_amps > 0, "no. of local amps");
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
	TEST_ASSERT(z == -1.9 + 0.8 * I, "amp=27, z=%f+%fi",
			creal(z), cimag(z));

	qreg_destroy(&reg);
}


void test_qreg_getsetamp_02(size_t tag)
{
	_Complex double z;
	struct qreg reg;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0,
		"cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z = rand_dbl() + rand_dbl() * I;
		AMPS[i] = z;
		qreg_setamp(&reg, i, z);
	}

	for (size_t i = 0; i < NUM_AMPS; i++) {
		qreg_getamp(&reg, i, &z);
		TEST_ASSERT(z == AMPS[i],
			"[%zu] i=%zu, z=%f+%fi, AMPS[i]=%f+%fi",
			tag, i, creal(z), cimag(z),
			creal(AMPS[i]), cimag(AMPS[i]));
	}

	qreg_destroy(&reg);
}

void test_qreg_zero(void)
{
	_Complex double z;
	struct qreg reg;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0,
		"cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z = rand_dbl() + rand_dbl() * I;
		AMPS[i] = z;
		qreg_setamp(&reg, i, z);
	}

	qreg_zero(&reg);

	for (size_t i = 0; i < NUM_AMPS; i++) {
		qreg_getamp(&reg, i, &z);
		TEST_ASSERT(z == 0.0,
			"i=%zu, z=%f+%fi", i, creal(z), cimag(z));
	}

	qreg_destroy(&reg);
}

/* 
 * Test if Pauli string equal to identity produces
 * just multiplication by phase.
 */
void test_qreg_paulirot_00(void)
{
	_Complex double z, z_exp;
	struct qreg reg;
	struct paulis ps_hi, ps_lo;
	double angle;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0,
		"cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z = rand_dbl() + rand_dbl() * I;
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
			"i=%zu, z=%f+%fi, z_exp=%f+%fi", i,
			creal(z), cimag(z), creal(z_exp), cimag(z_exp));
	}

	qreg_destroy(&reg);
}

/*
static void print_paulis(struct paulis ps)
{
	char buf[WIDTH+1];
	enum pauli_op op;
	for (size_t k = 0; k < WIDTH; k++) {
		op = paulis_get(ps, k);
		buf[k] = PAULI_LABEL[op];
	}
	buf[WIDTH] = '\0';
	printf("%s", buf);
}
*/

void test_qreg_paulirot_01(size_t tag)
{
	_Complex double z, z_exp;
	struct qreg reg;
	struct paulis ps, ps_hi, ps_lo;
	double angle;

	TEST_ASSERT(qreg_init(&reg, NUM_QUBITS) == 0,
		"cannot initialize qreg");
	for (size_t i = 0; i < NUM_AMPS; i++) {
		z = rand_dbl() + rand_dbl() * I;
		AMPS[i] = z;
		qreg_setamp(&reg, i, z);
	}

	ps = ps_hi = ps_lo = paulis_new();
	for (size_t k = 0; k < NUM_QUBITS; k++)
		paulis_set(&ps, rand_pauli_op(), k);
	paulis_split(ps, reg.qb_lo, reg.qb_hi, &ps_lo, &ps_hi);
	paulis_shr(&ps_hi, reg.qb_lo);
	angle = rand_dbl();
	qreg_paulirot(&reg, ps_hi, &ps_lo, &angle, 1);

	for (size_t i = 0; i < NUM_AMPS; i++) {
		_Complex double u = 1.0;
		
		size_t j = paulis_effect(ps, i, &u);
		z_exp = cos(angle) * AMPS[i] + I * u * sin(angle) * AMPS[j];
		
		qreg_getamp(&reg, i, &z);
		TEST_ASSERT(cabs(z - z_exp) < MARGIN,
			"[%zu] i=%zu, z=%f+%fi, z_exp=%f+%fi", tag, i,
			creal(z), cimag(z), creal(z_exp), cimag(z_exp));
	}

	qreg_destroy(&reg);
}

void TEST_MAIN(void)
{
	world_init((void *)0, (void *) 0);
	world_info(&WD);
	log_info("MPI world size: %d", WD.size);

	xoshiro256ss_init(&RNG, SEED);

	test_qreg_init();
	test_qreg_getsetamp_01();
	for (size_t k = 0; k < 9999; k++)
		test_qreg_getsetamp_02(k);
	test_qreg_zero();

	test_qreg_paulirot_00();
	for (size_t k = 0; k < 9999; k++)
		test_qreg_paulirot_01(k);

	world_fin();
}
