#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "algos/rayon.h"
#include "circ.h"
#include "data.h"

#define TEST(name, ...)                                                        \
	static int name(__VA_ARGS__)                                           \
	{                                                                      \
		const char *test_name = #name;                                 \
		int	    test_rc   = 0;

#define TEST_ASSERT(exp, ...)                                                  \
	if (!(exp)) {                                                          \
		fprintf(stderr, "FAILED %s: %s\n", test_name, #exp);           \
		fprintf(stderr, "  %s:%d\n  Reason: \"", __FILE__, __LINE__);  \
		fprintf(stderr, __VA_ARGS__);                                  \
		fprintf(stderr, "\"\n");                                       \
		goto test_error;                                               \
	}

#define TEST_FIN(STMNT)                                                        \
	goto test_exit;                                                        \
test_error:                                                                    \
	test_rc = -1;                                                          \
test_exit:                                                                     \
	STMNT;                                                                 \
	return test_rc;                                                        \
	}

#define TEST_CASE(exp) TEST_ASSERT(exp == 0, "%s", #exp)

#define MARGIN (0.0099)
static const char *CASE_DIR = DATA_DIR "/case-rand";

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

TEST(caserand, const char *prefix)
{
	static char buf[1024];
	snprintf(buf, 1024, "%s/%s", CASE_DIR, prefix);

	struct data dat, dat_ref;
	data_init(&dat);
	data_init(&dat_ref);
	TEST_ASSERT(read_data_ref(&dat, &dat_ref, buf) == 0, "Cannot read data "
							     "file");

	struct rayon_data rd;
	rayon_data_init(&rd);
	TEST_ASSERT(rayon_data_from_data(&rd, &dat) == 0, "Cannot parse "
							  "simulation data");
	TEST_ASSERT(rayon_simulate(&rd) == 0, "Simulation error");
	rayon_data_write_times(&dat.time_series, &rd.times);
	rayon_data_destroy(&rd);

	for (size_t i = 0; i < dat.time_series.num_steps; i++) {
		const _Complex double val = dat.time_series.values[i];
		const _Complex double ref = dat_ref.time_series.values[i];

		TEST_ASSERT(!(isnan(creal(val))),
			"Real value at index %zu is Nan", i);
		TEST_ASSERT(!(isnan(cimag(val))),
			"Imag value at index %zu is Nan", i);

		TEST_ASSERT(fabs(creal(val) - creal(ref)) < MARGIN,
			"Real diff exceeded margin (%f): val_re=%f, ref_re=%f",
			MARGIN, creal(val), creal(ref));
		TEST_ASSERT(fabs(cimag(val) - cimag(ref)) < MARGIN,
			"Imag diff exceeded margin (%f): val_re=%f, ref_re=%f",
			MARGIN, cimag(val), cimag(ref));
	}

	TEST_FIN({
		data_destroy(&dat);
		data_destroy(&dat_ref);
	});
}

TEST(caserand_suite, void)
{
	TEST_CASE(caserand("case-d9f603dc"));
	TEST_CASE(caserand("case-070d034c"));
	TEST_CASE(caserand("case-33427110"));
	TEST_CASE(caserand("case-28aa2595"));
	TEST_CASE(caserand("case-e1932ef1"));

	TEST_FIN();
}

int
main(void)
{
	return caserand_suite();
}