#include "c23_compat.h"
#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#include "phase2/circ.h"
#include "phase2/paulis.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0xb3a71f0e5c29d846)

static void t_hamil_init_free(void)
{
	struct circ_hamil hm;
	int rc;

	rc = circ_hamil_init(&hm, 5, 10);
	TEST_EQ(rc, 0);
	TEST_EQ(hm.qb, 5);
	TEST_EQ(hm.len, (size_t)10);
	TEST_ASSERT(hm.terms != nullptr, "terms should not be NULL");

	circ_hamil_free(&hm);
}

static void t_hamil_sort_lex(void)
{
	struct circ_hamil hm;
	int rc;

	rc = circ_hamil_init(&hm, 5, 5);
	TEST_EQ(rc, 0);

	/*
	 * Build 5 terms with distinct Pauli strings and known
	 * coefficients.  The coefficients serve as tags to verify
	 * that they move together with their operators during sort.
	 *
	 * Term 0: IIIZX  (qubit 0 = X, qubit 1 = Z)
	 * Term 1: IIIIX  (qubit 0 = X)
	 * Term 2: ZZIII  (qubit 3 = Z, qubit 4 = Z)
	 * Term 3: IIIIZ  (qubit 0 = Z)
	 * Term 4: IIIIY  (qubit 0 = Y)
	 */
	hm.terms[0].cf = 10.0;
	hm.terms[0].op = paulis_new();
	paulis_set(&hm.terms[0].op, PAULI_X, 0);
	paulis_set(&hm.terms[0].op, PAULI_Z, 1);

	hm.terms[1].cf = 20.0;
	hm.terms[1].op = paulis_new();
	paulis_set(&hm.terms[1].op, PAULI_X, 0);

	hm.terms[2].cf = 30.0;
	hm.terms[2].op = paulis_new();
	paulis_set(&hm.terms[2].op, PAULI_Z, 3);
	paulis_set(&hm.terms[2].op, PAULI_Z, 4);

	hm.terms[3].cf = 40.0;
	hm.terms[3].op = paulis_new();
	paulis_set(&hm.terms[3].op, PAULI_Z, 0);

	hm.terms[4].cf = 50.0;
	hm.terms[4].op = paulis_new();
	paulis_set(&hm.terms[4].op, PAULI_Y, 0);

	circ_hamil_sort_lex(&hm);

	/* Verify lexicographic order */
	for (size_t i = 0; i + 1 < hm.len; i++) {
		int cmp = paulis_cmp(hm.terms[i].op, hm.terms[i + 1].op);
		TEST_ASSERT(cmp <= 0,
			"terms[%zu] > terms[%zu] after sort (cmp=%d)",
			i, i + 1, cmp);
	}

	/*
	 * Verify coefficients moved with their operators.  Find
	 * each original operator and check its coefficient.
	 */
	for (size_t i = 0; i < hm.len; i++) {
		struct paulis op = hm.terms[i].op;
		double cf = hm.terms[i].cf;

		/* Reconstruct which original term this is by
		 * checking the Pauli string. */
		enum pauli_op p0 = paulis_get(op, 0);
		enum pauli_op p1 = paulis_get(op, 1);
		enum pauli_op p3 = paulis_get(op, 3);

		if (p0 == PAULI_X && p1 == PAULI_Z)
			TEST_EQ(cf, 10.0);
		else if (p0 == PAULI_X && p1 == PAULI_I && p3 == PAULI_I)
			TEST_EQ(cf, 20.0);
		else if (p0 == PAULI_I && p3 == PAULI_Z)
			TEST_EQ(cf, 30.0);
		else if (p0 == PAULI_Z)
			TEST_EQ(cf, 40.0);
		else if (p0 == PAULI_Y)
			TEST_EQ(cf, 50.0);
	}

	circ_hamil_free(&hm);
}

static void t_muldet_init_free(void)
{
	struct circ_muldet md;
	int rc;

	rc = circ_muldet_init(&md, 3);
	TEST_EQ(rc, 0);
	TEST_EQ(md.len, (size_t)3);
	TEST_ASSERT(md.dets != nullptr, "dets should not be NULL");

	circ_muldet_free(&md);
}

static void t_values_init_free(void)
{
	struct circ_values vals;
	int rc;

	rc = circ_values_init(&vals, 7);
	TEST_EQ(rc, 0);
	TEST_EQ(vals.len, (size_t)7);
	TEST_ASSERT(vals.z != nullptr, "z should not be NULL");

	/* Write known complex values, then read back and verify. */
	for (size_t k = 0; k < 7; k++)
		vals.z[k] = (double)k + (double)(k + 1) * I;

	for (size_t k = 0; k < 7; k++) {
		_Complex double expected = (double)k + (double)(k + 1) * I;
		TEST_ASSERT(vals.z[k] == expected,
			"vals.z[%zu] mismatch", k);
	}

	circ_values_free(&vals);
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	t_hamil_init_free();
	t_hamil_sort_lex();
	t_muldet_init_free();
	t_values_init_free();

	world_free();
}
