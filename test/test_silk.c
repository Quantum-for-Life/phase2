#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "algos/silk.h"
#include "data2.h"

#include "test.h"

#define MARGIN (0.0099)
static const char *CASE_DIR = PH2_SIMUL_DATA "/case-rand";

static int caserand(const char *prefix)
{
	char filename[1024] = { 0 };

	snprintf(filename, 1024, "%s/%s.h5", CASE_DIR, prefix);
	data2_id fid = data2_open(filename);
	if (fid == DATA2_INVALID_FID) {
		TEST_FAIL("Cannot read data file: %s", filename);
		goto err_data_open;
	}

	struct silk_data rd;
	silk_data_init(&rd);
	if (silk_data_from_data(&rd, fid) != 0) {
		TEST_FAIL("Cannot parse simulation data");
		goto err_rd_read;
	}

	if (silk_simulate(&rd) != 0) {
		TEST_FAIL("Simulation error");
		goto err_simul;
	}

	snprintf(filename, 1024, "%s/%s.h5_solved", CASE_DIR, prefix);
	data2_id fid_ref = data2_open(filename);
	if (fid_ref == DATA2_INVALID_FID) {
		TEST_FAIL("Cannot read data file: %s", filename);
		goto err_data_open_ref;
	}

	struct silk_data rd_ref;
	silk_data_init(&rd_ref);
	if (silk_data_from_data(&rd_ref, fid_ref) != 0) {
		TEST_FAIL("Cannot parse reference data");
		goto err_rd_ref_read;
	}

	for (size_t i = 0; i < rd.times.num_steps; i++) {
		const _Complex double val = rd.times.steps[i].val;
		const _Complex double ref = rd_ref.times.steps[i].val;

		if (isnan(creal(val))) {
			TEST_FAIL("Real value at index %zu is Nan", i);
			goto error;
		}
		if (isnan(cimag(val))) {
			TEST_FAIL("Imag value at index %zu is Nan", i);
			goto error;
		}
		if (fabs(creal(val) - creal(ref)) > MARGIN) {
			TEST_FAIL("Real diff exceeded margin (%f): val_re=%f, "
				  "ref_re=%f",
				MARGIN, creal(val), creal(ref));
			goto error;
		}
		if (fabs(cimag(val) - cimag(ref)) > MARGIN) {
			TEST_FAIL(
				"Imag diff exceeded margin (%f): val_re=%f, ref_re=%f",
				MARGIN, cimag(val), cimag(ref));
			goto error;
		}
	}

	silk_data_destroy(&rd_ref);
	data2_close(fid_ref);
	silk_data_destroy(&rd);
	data2_close(fid);
	return 0;

error:
	silk_data_destroy(&rd_ref);
err_rd_ref_read:
	data2_close(fid_ref);
err_data_open_ref:
err_simul:
	silk_data_destroy(&rd);
err_rd_read:
	data2_close(fid);
err_data_open:
	return -1;
}

int main(void)
{
	if (caserand("case-d9f603dc") != 0)
		goto error;
	if (caserand("case-070d034c") != 0)
		goto error;
	if (caserand("case-33427110") != 0)
		goto error;
	if (caserand("case-28aa2595") != 0)
		goto error;
	if (caserand("case-e1932ef1") != 0)
		goto error;

	return 0;
error:
	return -1;
}
