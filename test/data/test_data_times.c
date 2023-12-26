#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "data.h"

#include "test.h"
#include "test_data.h"

#define MARGIN (10e-5)

static int
test_getnums(void)
{
	int rc = 0;
	for (size_t i = 0; i < NUM_TEST_FILES; i++) {
		struct test_data td	  = TEST_DATA[i];
		const char	*filename = td.filename;

		data_id fid = data2_open(filename);
		if (fid == DATA_INVALID_FID) {
			TEST_FAIL("open file: %s", filename);
			rc = -1;
			break;
		}

		size_t num_steps = 0;
		if (data2_times_getnums(fid, &num_steps) < 0) {
			TEST_FAIL("read times getnums()");
			rc = -1;
			break;
		}
		if (num_steps != td.num_steps) {
			TEST_FAIL("wrong number of steps in time series: %zu",
				num_steps);
			rc = -1;
		}
		data2_close(fid);
	}

	return rc;
}

static int
iter_count(double t, _Complex double v, void *op_data)
{
	(void)t;
	(void)v;

	int *count = op_data;
	(*count)++;

	return 0;
}

static int
iter_count_onlytwo(double t, _Complex double v, void *op_data)
{
	(void)t;
	(void)v;

	int *count = op_data;
	(*count)++;
	if (*count == 2)
		return 77;

	return 0;
}

struct iter_store_data {
	size_t		idx;
	size_t		num_steps;
	double		times[128];
	_Complex double values[128];
};

static int
iter_store(double t, _Complex double v, void *op_data)
{
	struct iter_store_data *idat = op_data;

	idat->times[idat->idx]	= t;
	idat->values[idat->idx] = v;
	idat->idx++;

	return 0;
}

static int
test_iter0(void)
{
	const struct test_data td  = TEST_DATA[0];
	data_id		       fid = data2_open(td.filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("open file: %s", td.filename);
		return -1;
	}

	size_t num_steps;
	if (data2_times_getnums(fid, &num_steps) < 0) {
		TEST_FAIL("wron number of steps: %zu", num_steps);
		goto err;
	}

	int count = 0;
	if (data2_times_foreach(fid, iter_count, &count) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	if (count != td.num_steps) {
		TEST_FAIL("number of iterations: %d", count);
		goto err;
	}

	count = 0;
	if (data2_times_foreach(fid, iter_count_onlytwo, &count) == 0) {
		TEST_FAIL("iteration didn't terminate early");
		goto err;
	}
	if (count != 2) {
		TEST_FAIL("number of iterations: %d", count);
		goto err;
	}

	struct iter_store_data idat = { .idx = 0, .num_steps = num_steps };
	if (data2_times_foreach(fid, iter_store, &idat) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	double		exp_times[4]  = { 206, 125, 138, 21 };
	_Complex double exp_values[4] = {
		NAN + _Complex_I * NAN,
		NAN + _Complex_I * NAN,
		NAN + _Complex_I * NAN,
		NAN + _Complex_I * NAN,
	};
	for (size_t i = 0; i < num_steps; i++) {
		if (fabs(idat.times[i] - exp_times[i]) > MARGIN) {
			TEST_FAIL("time %f vs. expected: %f", idat.times[i],
				exp_times[i]);
			goto err;
		}
		if (cabs(idat.values[i] - exp_values[i]) > MARGIN) {
			TEST_FAIL("value %f + %f i vs. expected: %f + %f i",
				creal(idat.values[i]), cimag(idat.values[i]),
				creal(exp_values[i]), cimag(exp_values[i]));
			goto err;
		}
	}

	data2_close(fid);
	return 0;
err:
	data2_close(fid);
	return -1;
}

static int
test_iter1(void)
{
	const struct test_data td  = TEST_DATA[1];
	data_id		       fid = data2_open(td.filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("open file: %s", td.filename);
		return -1;
	}

	size_t num_steps;
	if (data2_times_getnums(fid, &num_steps) < 0) {
		TEST_FAIL("wron number of steps: %zu", num_steps);
		goto err;
	}

	int count = 0;
	if (data2_times_foreach(fid, iter_count, &count) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	if (count != td.num_steps) {
		TEST_FAIL("number of iterations: %d", count);
		goto err;
	}

	count = 0;
	if (data2_times_foreach(fid, iter_count_onlytwo, &count) == 0) {
		TEST_FAIL("iteration didn't terminate early");
		goto err;
	}
	if (count != 2) {
		TEST_FAIL("number of iterations: %d", count);
		goto err;
	}

	struct iter_store_data idat = { .idx = 0, .num_steps = num_steps };
	if (data2_times_foreach(fid, iter_store, &idat) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	double exp_times[16] = { 20, 20, 21, 10, 11, 37, 48, 9, 57, 41, 20, 30,
		21, 46, 31, 52 };
	_Complex double exp_values[16] = { 0.034389 + _Complex_I * 0.209985,
		0.034389 + _Complex_I * 0.209985,
		0.106208 + _Complex_I * 0.188456,
		0.00669263 + _Complex_I * -0.187835,
		0.141371 + _Complex_I * -0.147162,
		-0.561908 + _Complex_I * 0.0917758,
		0.295048 + _Complex_I * 0.0561566,
		-0.157683 + _Complex_I * -0.195077,
		0.0543732 + _Complex_I * -0.242585,
		-0.714111 + _Complex_I * -0.12588,
		0.034389 + _Complex_I * 0.209985,
		0.0471345 + _Complex_I * -0.370453,
		0.106208 + _Complex_I * 0.188456,
		0.472685 + _Complex_I * 0.0128869,
		0.198564 + _Complex_I * -0.319754,
		-0.25097 + _Complex_I * -0.146826 };
	for (size_t i = 0; i < num_steps; i++) {
		if (fabs(idat.times[i] - exp_times[i]) > MARGIN) {
			TEST_FAIL("time %f vs. expected: %f", idat.times[i],
				exp_times[i]);
			goto err;
		}
		if (cabs(idat.values[i] - exp_values[i]) > MARGIN) {
			TEST_FAIL("value %f + %f i vs. expected: %f + %f i",
				creal(idat.values[i]), cimag(idat.values[i]),
				creal(exp_values[i]), cimag(exp_values[i]));
			goto err;
		}
	}

	data2_close(fid);
	return 0;
err:
	data2_close(fid);
	return -1;
}

int
test_data_times()
{
	if (test_getnums() < 0)
		goto err;
	if (test_iter0() < 0)
		goto err;
	if (test_iter1() < 0)
		goto err;

	return 0;
err:
	return -1;
}
