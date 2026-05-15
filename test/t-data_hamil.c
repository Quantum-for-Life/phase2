/*
 * Test the circ_hamil_load API on every committed fixture
 * carrying a /pauli_hamil group.  Asserts: qb, len, and (for
 * the second fixture) per-term coefficient and packed Pauli
 * operator.  Coefficients are scaled by the on-disk norm at
 * load time; the per-term cf check folds the norm in
 * implicitly.
 */
#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "phase2/circ.h"
#include "ph2run/data.h"
#include "phase2/paulis.h"
#include "phase2/world.h"

#include "t-data.h"
#include "test.h"

#define WD_SEED UINT64_C(0xc9c70166d249f2d4)

#define MARGIN (1.0e-8)

static int t_dims_and_packing(void)
{
	int rc = 0;
	for (size_t i = 0; i < NUM_TEST_FILES; i++) {
		const struct test_data td = TEST_DATA[i];

		data_id fid = data_open(td.filename);
		if (fid == DATA_INVALID_FID) {
			TEST_FAIL("open file: %s", td.filename);
			rc = -1;
			break;
		}

		struct circ_hamil hm;
		if (circ_hamil_load(fid, &hm) < 0) {
			TEST_FAIL("circ_hamil_load on %s", td.filename);
			rc = -1;
			data_close(fid);
			break;
		}
		if (hm.qb != td.num_qubits) {
			TEST_FAIL("wrong number of qubits: %zu vs %zu",
				(size_t)hm.qb, td.num_qubits);
			rc = -1;
		}
		if (hm.len != td.num_terms) {
			TEST_FAIL("wrong number of terms: %zu vs %zu",
				hm.len, td.num_terms);
			rc = -1;
		}
		circ_hamil_free(&hm);
		data_close(fid);
	}
	return rc;
}

static int t_terms(void)
{
	const struct test_data td = TEST_DATA[1];
	data_id fid = data_open(td.filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("open file: %s", td.filename);
		return -1;
	}

	struct circ_hamil hm;
	if (circ_hamil_load(fid, &hm) < 0) {
		TEST_FAIL("circ_hamil_load on %s", td.filename);
		data_close(fid);
		return -1;
	}

	int rc = 0;
	const double exp_coeff[] = { 0.64604963, 0.16592673, 0.90327525,
		-0.18683327, -0.18315831, 0.57830137, 0.71210119, -0.96550733,
		0.21017606, -0.84378561 };
	const unsigned char exp_paulis[] = { 0, 3, 1, 3, 0, 3, 0, 1, 0, 1, 1,
		3, 0, 3, 0, 2, 1, 3, 1, 0, 0, 3, 2, 1, 2, 3, 1, 1, 1, 1 };
	for (size_t i = 0; i < hm.len; i++) {
		const double exp_cf = exp_coeff[i] * td.norm;
		if (fabs(hm.terms[i].cf - exp_cf) > MARGIN) {
			TEST_FAIL("term[%zu].cf: %f vs %f", i,
				hm.terms[i].cf, exp_cf);
			rc = -1;
		}
		struct paulis exp_op = paulis_new();
		for (uint32_t j = 0; j < hm.qb; j++)
			paulis_set(&exp_op,
				exp_paulis[i * hm.qb + j], j);
		if (!paulis_eq(hm.terms[i].op, exp_op)) {
			TEST_FAIL("term[%zu].op mismatch", i);
			rc = -1;
		}
	}
	circ_hamil_free(&hm);
	data_close(fid);
	return rc;
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	if (t_dims_and_packing() < 0)
		TEST_FAIL("dims_and_packing");
	if (t_terms() < 0)
		TEST_FAIL("terms");

	world_free();
}
