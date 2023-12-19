#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "algos/rayon.h"
#include "circ.h"
#include "data.h"

#include "test.h"

#define MARGIN (0.0099)
static const char *CASE_DIR = PH2_SIMUL_DATA "/case-rand";

int
datio_read_file(struct data *dat, const char *filename)
{
	const data_id fid = data_file_open(filename);
	const int     rc  = data_parse(dat, fid);
	data_file_close(fid);

	return rc;
}

int
read_data_ref(struct data *dat, struct data *dat_ref, const char *path_prefix)
{
	const size_t buf_len = strlen(path_prefix) + strlen(".h5_solved") + 1;
	char	    *buf     = calloc(buf_len, sizeof(*buf));
	if (!buf)
		return -1;

	snprintf(buf, buf_len, "%s.h5", path_prefix);
	if (datio_read_file(dat, buf) < 0)
		goto error;
	snprintf(buf, buf_len, "%s.h5_solved", path_prefix);
	if (datio_read_file(dat_ref, buf) < 0)
		goto error;

	free(buf);
	return 0;
error:
	free(buf);
	return -1;
}

static int
caserand(const char *prefix)
{
	static char buf[1024];
	snprintf(buf, 1024, "%s/%s", CASE_DIR, prefix);

	struct data dat, dat_ref;
	data_init(&dat);
	data_init(&dat_ref);
	if (read_data_ref(&dat, &dat_ref, buf) != 0) {
		TEST_FAIL("Cannot read data file");
		goto error;
	}

	struct rayon_data rd;
	rayon_data_init(&rd);
	if (rayon_data_from_data(&rd, &dat) != 0) {
		TEST_FAIL("Cannot parse simulation data");
		goto error;
	}
	if (rayon_simulate(&rd) != 0) {
		TEST_FAIL("Simulation error");
		goto error;
	}
	rayon_data_write_times(&dat.time_series, &rd.times);
	rayon_data_destroy(&rd);

	for (size_t i = 0; i < dat.time_series.num_steps; i++) {
		const _Complex double val = dat.time_series.values[i];
		const _Complex double ref = dat_ref.time_series.values[i];

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

	data_destroy(&dat);
	data_destroy(&dat_ref);
	return 0;
error:
	data_destroy(&dat);
	data_destroy(&dat_ref);
	return -1;
}

int
main(void)
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
