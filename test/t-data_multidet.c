/*
 * Test circ_muldet_load on every committed fixture carrying
 * a /state_prep/multidet group.  Asserts: packed length, that
 * each packed idx fits in the fixture's qubit count, and (for
 * the second fixture) the exact contents of the basis-state
 * indices and complex coefficients.
 */
#include "c23_compat.h"
#include <complex.h>
#include <stdint.h>
#include <stdio.h>

#include "phase2/circ.h"
#include "ph2run/data.h"
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

		struct circ_muldet md;
		if (circ_muldet_load(fid, &md) < 0) {
			TEST_FAIL("circ_muldet_load on %s", td.filename);
			rc = -1;
			data_close(fid);
			break;
		}
		if (md.len != td.num_dets) {
			TEST_FAIL("wrong number of dets: %zu vs %zu",
				md.len, td.num_dets);
			rc = -1;
		}
		const uint64_t qb_mask = (td.num_qubits >= 64)
			? UINT64_MAX
			: ((UINT64_C(1) << td.num_qubits) - 1);
		for (size_t j = 0; j < md.len; j++) {
			if (md.dets[j].idx & ~qb_mask) {
				TEST_FAIL("idx[%zu]=0x%llx exceeds %zu qubits",
					j,
					(unsigned long long)md.dets[j].idx,
					td.num_qubits);
				rc = -1;
			}
		}
		circ_muldet_free(&md);
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

	struct circ_muldet md;
	if (circ_muldet_load(fid, &md) < 0) {
		TEST_FAIL("circ_muldet_load");
		data_close(fid);
		return -1;
	}

	int rc = 0;
	const uint64_t exp_idx[] = { 4, 5, 6 };
	const _Complex double exp_coeff[] = { CMPLX(0.108292, 0.333811),
		CMPLX(0.0491404, 0.613936), CMPLX(0.565802, 0.421163) };

	if (md.len != 3) {
		TEST_FAIL("expected 3 dets, got %zu", md.len);
		rc = -1;
	} else {
		for (size_t i = 0; i < md.len; i++) {
			if (md.dets[i].idx != exp_idx[i]) {
				TEST_FAIL("idx[%zu]: %lu vs %lu",
					i, (unsigned long)md.dets[i].idx,
					(unsigned long)exp_idx[i]);
				rc = -1;
			}
			const _Complex double cf = md.dets[i].cf;
			if (cabs(cf - exp_coeff[i]) > MARGIN) {
				TEST_FAIL("coeff[%zu]: %f+%fi vs %f+%fi",
					i, creal(cf), cimag(cf),
					creal(exp_coeff[i]),
					cimag(exp_coeff[i]));
				rc = -1;
			}
		}
	}

	circ_muldet_free(&md);
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
