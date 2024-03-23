#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ev.h"
#include "paulis.h"
#include "qreg.h"
#include "types.h"

#define NUM_QUBITS (20)
#define NUM_AMPS ((u64)1 << NUM_QUBITS)

static int RET = 0;

#define TEST_FAIL(...)                                                         \
	{                                                                      \
		fprintf(stderr, "FAILED %s:%d \"", __FILE__, __LINE__);        \
		fprintf(stderr, __VA_ARGS__);                                  \
		fprintf(stderr, "\"\n");                                       \
		RET = -1;                                                      \
	}

static struct ev   EV;
static struct qreg QREG;

void test_rot0(void)
{
	qreg_zerostate(&QREG);

	const struct paulis code_hi = paulis_new();
	const struct paulis code_lo = paulis_new();
	const f64	    angle   = M_PI_2;
	qreg_paulirot(&QREG, code_hi, &code_lo, &angle, 1);

	f64 amp[2];
	qreg_getamp(&QREG, 0, &amp);
	if (abs(amp[0]) > DBL_EPSILON || abs(amp[1] - 1.0) > DBL_EPSILON)
		TEST_FAIL("amp[0] should be 0.0+1.0i");

	for (u64 i = 1; i < NUM_AMPS; i++) {
		qreg_getamp(&QREG, i, &amp);
		if (abs(amp[0]) > DBL_EPSILON || abs(amp[1] > DBL_EPSILON)) {
			TEST_FAIL("amp[%lu] should be 0.0", i);
			break;
		}
	}
}

void test_rot1(void)
{
	qreg_zerostate(&QREG);

	struct paulis code = paulis_new();
	for (u32 i = 0; i < NUM_QUBITS; i++)
		paulis_set(&code, PAULI_X, i);

	struct paulis code_lo, code_hi;
	paulis_split(code, QREG.qb_lo, QREG.qb_hi, &code_lo, &code_hi);
	paulis_shr(&code_hi, QREG.qb_lo);

	const f64 angle = M_PI_2;
	qreg_paulirot(&QREG, code_hi, &code_lo, &angle, 1);

	f64 amp[2];
	qreg_getamp(&QREG, 0, &amp);
	for (u64 i = 0; i < NUM_AMPS - 1; i++) {
		qreg_getamp(&QREG, i, &amp);
		if (abs(amp[0]) > DBL_EPSILON || abs(amp[1] > DBL_EPSILON)) {
			TEST_FAIL("amp[%lu] should be 0.0, but is %f+%fi", i,
				amp[0], amp[1]);
			break;
		}
	}

	qreg_getamp(&QREG, NUM_AMPS - 1, &amp);
	if (abs(amp[0]) > DBL_EPSILON || abs(amp[1] - 1.0) > DBL_EPSILON)
		TEST_FAIL("amp[%lu] should be 0.0+1.0i, but is %f+%fi",
			NUM_AMPS - 1, amp[0], amp[1]);
}

void cleanup(void)
{
	qreg_destroy(&QREG);
	ev_destroy(&EV);
}

int main(void)
{
	if (ev_init(&EV) < 0) TEST_FAIL("init environment");
	if (qreg_init(&QREG, NUM_QUBITS, &EV) < 0) TEST_FAIL("init qreg");
	atexit(cleanup);

	test_rot0();
	test_rot1();

	return RET;
}
