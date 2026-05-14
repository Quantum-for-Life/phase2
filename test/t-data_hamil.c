/*
 * Test the data_hamil_load API on every committed fixture
 * carrying a /pauli_hamil group.  Asserts: nqb, nterms,
 * norm, and (for the second fixture) the exact contents of
 * the coefficient + Pauli arrays.
 */
#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "phase2/data.h"
#include "phase2/world.h"

#include "t-data.h"
#include "test.h"

#define WD_SEED UINT64_C(0xc9c70166d249f2d4)

#define MARGIN (1.0e-8)

static int t_dims_and_norm(void)
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

		struct data_hamil h;
		if (data_hamil_load(fid, &h) < 0) {
			TEST_FAIL("data_hamil_load on %s", td.filename);
			rc = -1;
			data_close(fid);
			break;
		}
		if (h.nqb != td.num_qubits) {
			TEST_FAIL("wrong number of qubits: %zu vs %zu",
				(size_t)h.nqb, td.num_qubits);
			rc = -1;
		}
		if (h.nterms != td.num_terms) {
			TEST_FAIL("wrong number of terms: %zu vs %zu",
				h.nterms, td.num_terms);
			rc = -1;
		}
		if (fabs(h.norm - td.norm) > MARGIN) {
			TEST_FAIL("wrong norm: %f vs %f", h.norm, td.norm);
			rc = -1;
		}
		data_hamil_free(&h);
		data_close(fid);
	}
	return rc;
}

static int t_arrays(void)
{
	const struct test_data td = TEST_DATA[1];
	data_id fid = data_open(td.filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("open file: %s", td.filename);
		return -1;
	}

	struct data_hamil h;
	if (data_hamil_load(fid, &h) < 0) {
		TEST_FAIL("data_hamil_load on %s", td.filename);
		data_close(fid);
		return -1;
	}

	int rc = 0;
	const double exp_coeff[] = { 0.64604963, 0.16592673, 0.90327525,
		-0.18683327, -0.18315831, 0.57830137, 0.71210119, -0.96550733,
		0.21017606, -0.84378561 };
	const unsigned char exp_paulis[] = { 0, 3, 1, 3, 0, 3, 0, 1, 0, 1, 1,
		3, 0, 3, 0, 2, 1, 3, 1, 0, 0, 3, 2, 1, 2, 3, 1, 1, 1, 1 };
	for (size_t i = 0; i < h.nterms; i++) {
		if (fabs(h.cfs[i] - exp_coeff[i]) > MARGIN) {
			TEST_FAIL("coeff[%zu]: %f vs %f", i, h.cfs[i],
				exp_coeff[i]);
			rc = -1;
		}
	}
	for (size_t i = 0; i < h.nterms * h.nqb; i++) {
		if (h.paulis[i] != exp_paulis[i]) {
			TEST_FAIL("pauli[%zu]: %u vs %u", i, h.paulis[i],
				exp_paulis[i]);
			rc = -1;
		}
	}
	data_hamil_free(&h);
	data_close(fid);
	return rc;
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	if (t_dims_and_norm() < 0)
		TEST_FAIL("dims_and_norm");
	if (t_arrays() < 0)
		TEST_FAIL("arrays");

	world_free();
}
