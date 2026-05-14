/*
 * Test data_multidet_load on every committed fixture
 * carrying a /state_prep/multidet group.  Asserts: nqb,
 * ndets, and (for the second fixture) the exact contents
 * of the basis-state indices and complex coefficients.
 */
#include "c23_compat.h"
#include <complex.h>
#include <stdint.h>
#include <stdio.h>

#include "phase2/data.h"
#include "phase2/world.h"

#include "t-data.h"
#include "test.h"

#define WD_SEED UINT64_C(0x8adaececa772f40d)

#define MARGIN (1.0e-6)

static int t_dims(void)
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

		struct data_multidet m;
		if (data_multidet_load(fid, &m) < 0) {
			TEST_FAIL("data_multidet_load on %s", td.filename);
			rc = -1;
			data_close(fid);
			break;
		}
		if (m.nqb != td.num_qubits) {
			TEST_FAIL("wrong number of qubits: %zu vs %zu",
				(size_t)m.nqb, td.num_qubits);
			rc = -1;
		}
		if (m.ndets != td.num_dets) {
			TEST_FAIL("wrong number of dets: %zu vs %zu",
				m.ndets, td.num_dets);
			rc = -1;
		}
		data_multidet_free(&m);
		data_close(fid);
	}
	return rc;
}

/* Build the basis-state index (LSB = qubit 0) from a row of
 * the raw `dets` byte array, the way circ.c unpacks it. */
static uint64_t pack_idx(const unsigned char *row, uint32_t nqb)
{
	uint64_t idx = 0;
	for (uint32_t j = 0; j < nqb; j++)
		idx += (uint64_t)row[j] << j;
	return idx;
}

static int t_arrays(void)
{
	const struct test_data td = TEST_DATA[1];
	data_id fid = data_open(td.filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("open file: %s", td.filename);
		return -1;
	}

	struct data_multidet m;
	if (data_multidet_load(fid, &m) < 0) {
		TEST_FAIL("data_multidet_load");
		data_close(fid);
		return -1;
	}

	int rc = 0;
	const uint64_t exp_idx[] = { 4, 5, 6 };
	const _Complex double exp_coeff[] = { CMPLX(0.108292, 0.333811),
		CMPLX(0.0491404, 0.613936), CMPLX(0.565802, 0.421163) };

	if (m.ndets != 3) {
		TEST_FAIL("expected 3 dets, got %zu", m.ndets);
		rc = -1;
	} else {
		for (size_t i = 0; i < m.ndets; i++) {
			const uint64_t idx = pack_idx(
				&m.dets[i * m.nqb], m.nqb);
			if (idx != exp_idx[i]) {
				TEST_FAIL("idx[%zu]: %lu vs %lu",
					i, (unsigned long)idx,
					(unsigned long)exp_idx[i]);
				rc = -1;
			}
			const _Complex double cf = CMPLX(
				m.cfs[2 * i], m.cfs[2 * i + 1]);
			if (cabs(cf - exp_coeff[i]) > MARGIN) {
				TEST_FAIL("coeff[%zu]: %f+%fi vs %f+%fi",
					i, creal(cf), cimag(cf),
					creal(exp_coeff[i]),
					cimag(exp_coeff[i]));
				rc = -1;
			}
		}
	}

	data_multidet_free(&m);
	data_close(fid);
	return rc;
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	if (t_dims() < 0)
		TEST_FAIL("dims");
	if (t_arrays() < 0)
		TEST_FAIL("arrays");

	world_free();
}
