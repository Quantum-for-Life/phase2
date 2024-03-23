#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "algos/silk.h"
#include "data2.h"

#include "test.h"

#define NUM_STEPS (128)

#define NUM_CASES (16)
const char *CASEIDS[NUM_CASES] = {

	"19f7187e", "f2a421c7", "eedff5bd", "04ba000b", "e502d0bb", "caf81b88",
	"beb0578e", "39663512", "00358789", "cf50c04b", "6171a3c8", "fb20d1f9",
	"39a48d76", "2779488a", "e13d2e13", "019c97dc"

};

#define MARGIN (1e-14)
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
	silk_data_init(&rd, NUM_STEPS);
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

	_Complex double ref_values[NUM_STEPS];
	if (data2_trotter_read_values_test(fid_ref, ref_values) != 0) {
		TEST_FAIL("Cannot parse reference data");
		goto err_rd_ref_read;
	}

	for (size_t i = 0; i < rd.num_steps; i++) {
		const _Complex double val = rd.trotter_steps[i];
		const _Complex double ref = ref_values[i];

		if (fabs(creal(val) - creal(ref)) > MARGIN) {
			TEST_FAIL("Real diff exceeded margin (%.16f): "
				  "val_re=%.16f, ref_re=%.16f",
				MARGIN, creal(val), creal(ref));
			goto error;
		}
		if (fabs(cimag(val) - cimag(ref)) > MARGIN) {
			TEST_FAIL("Imag diff exceeded margin (%.16f): "
				  "val_re=%.16f, ref_re=%.16f",
				MARGIN, cimag(val), cimag(ref));
			goto error;
		}
	}

	data2_close(fid_ref);
	silk_data_destroy(&rd);
	data2_close(fid);
	return 0;

error:
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
	for (size_t i = 0; i < NUM_CASES; i++) {
		char buf[64];
		snprintf(buf, 64, "case-%s", CASEIDS[i]);
		if (caserand(buf) != 0) {
			TEST_FAIL("random case %s", buf);
			goto error;
		}
	}

	return 0;
error:
	return -1;
}
