#include <complex.h>
#include <stdio.h>

#include "phase2/data.h"
#include "phase2/world.h"

#include "t-data.h"
#include "test.h"

#define MARGIN (1.0e-6)

static int t_get_nums(void)
{
	int rc = 0;
	for (size_t i = 0; i < NUM_TEST_FILES; i++) {
		struct test_data td = TEST_DATA[i];
		const char *filename = td.filename;

		data_id fid = data_open(filename);
		if (fid == DATA_INVALID_FID) {
			TEST_FAIL("open file: %s", filename);
			rc = -1;
			break;
		}

		size_t num_qubits = 0, num_dets = 0;
		if (data_multidet_getnums(fid, &num_qubits, &num_dets) < 0) {
			TEST_FAIL("read multidet getnums()");
			rc = -1;
			break;
		}
		if (num_qubits != td.num_qubits) {
			TEST_FAIL("wrong number of qubits: %zu", num_qubits);
			rc = -1;
		}
		if (num_dets != td.num_dets) {
			TEST_FAIL("wrong number of dets: %zu", num_dets);
			rc = -1;
		}
		data_close(fid);
	}

	return rc;
}

static int iter_count_dets(double coeff[2], size_t idx, void *op_data)
{
	(void)coeff;
	(void)idx;

	int *count = op_data;
	(*count)++;

	return 0;
}

static int iter_count_dets_onlytwo(double coeff[2], size_t idx, void *op_data)
{
	(void)coeff;
	(void)idx;

	int *count = op_data;
	(*count)++;

	if (*count == 2)
		return 91;

	return 0;
}

struct iter_store {
	size_t index;
	_Complex double coeff[128];
	size_t idx[128];
};

static int iter_store_dets(double coeff[2], const size_t idx, void *op_data)
{
	struct iter_store *is = op_data;

	is->coeff[is->index] = coeff[0] + _Complex_I * coeff[1];
	is->idx[is->index] = idx;
	is->index++;

	return 0;
}

static int t_iter0(void)
{
	const struct test_data td = TEST_DATA[0];
	data_id fid = data_open(td.filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("open file: %s", td.filename);
		return -1;
	}

	size_t num_qubits, num_dets;
	data_multidet_getnums(fid, &num_qubits, &num_dets);

	int count = 0;
	if (data_multidet_foreach(fid, iter_count_dets, &count) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	if (count != 1) {
		TEST_FAIL("number of iterations: %d", count);
		goto err;
	}
	count = 0;
	if (data_multidet_foreach(fid, iter_count_dets_onlytwo, &count) != 0) {
		TEST_FAIL("iteration return wrong code");
		goto err;
	}
	if (count != 1) {
		TEST_FAIL("number of iterations: %d", count);
		goto err;
	}

	data_close(fid);
	return 0;
err:
	data_close(fid);
	return -1;
}

static int t_iter1(void)
{
	const struct test_data td = TEST_DATA[1];
	data_id fid = data_open(td.filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("open file: %s", td.filename);
		return -1;
	}

	size_t num_qubits, num_dets;
	data_multidet_getnums(fid, &num_qubits, &num_dets);

	int count = 0;
	if (data_multidet_foreach(fid, iter_count_dets, &count) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	if (count != 3) {
		TEST_FAIL("number of iterations: %d", count);
		goto err;
	}
	count = 0;
	if (data_multidet_foreach(fid, iter_count_dets_onlytwo, &count) != 91) {
		TEST_FAIL("iteration returned wrong code");
		goto err;
	}
	if (count != 2) {
		TEST_FAIL("number of iterations: %d", count);
		goto err;
	}

	struct iter_store is;
	size_t exp_idx[] = { 4, 5, 6 };
	_Complex double exp_coeff[] = { 0.108292 + I * 0.333811,
		0.0491404 + I * 0.613936, 0.565802 + I * 0.421163 };

	is.index = 0;
	if (data_multidet_foreach(fid, iter_store_dets, &is) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	if (is.index != 3) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}

	for (size_t i = 0; i < td.num_dets; i++) {
		if (exp_idx[i] != is.idx[i]) {
			TEST_FAIL("idx cmp fail: %zu vs. %zu", exp_idx[i],
				is.idx[i]);
			goto err;
		}
		if (cabs(exp_coeff[i] - is.coeff[i]) > MARGIN) {
			TEST_FAIL("complex cmp fail %f + %f i vs. %f + %f i",
				creal(is.coeff[i]), cimag(is.coeff[i]),
				creal(exp_coeff[i]), cimag(exp_coeff[i]));
			goto err;
		}
	}

	data_close(fid);
	return 0;
err:
	data_close(fid);
	return -1;
}

static void TEST_MAIN(void)
{
	world_init((void *)0, (void *) 0);

	if (t_get_nums() < 0) {
		TEST_FAIL("getnums");
		return;
	}
	if (t_iter0() < 0) {
		TEST_FAIL("iter0");
		return;
	}
	if (t_iter1() < 0) {
		TEST_FAIL("iter1");
		return;
	}

	world_fin();
}
