#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "data.h"

#include "test.h"
#include "test-data.h"

#define MARGIN (10e-8)

static int test_getnums(void)
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

		size_t num_qubits = 0, num_terms = 0;
		if (data_hamil_getnums(fid, &num_qubits, &num_terms) < 0) {
			TEST_FAIL("read hamil getnums()");
			rc = -1;
			break;
		}
		if (num_qubits != td.num_qubits) {
			TEST_FAIL("wrong number of qubits in hamil: %zu",
				num_qubits);
			rc = -1;
		}
		if (num_terms != td.num_terms) {
			TEST_FAIL(
				"wrong number of hamil terms: %zu", num_terms);
			rc = -1;
		}
		data_close(fid);
	}

	return rc;
}

static int test_getnorm(void)
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

		double norm;
		if (data_hamil_getnorm(fid, &norm) < 0) {
			TEST_FAIL("read multidet getnums()");
			rc = -1;
			break;
		}
		if (fabs(norm - td.norm) > MARGIN) {
			TEST_FAIL("norm for hamil: %f (expect. %f)", norm,
				td.norm);
			rc = -1;
		}

		data_close(fid);
	}

	return rc;
}

static int iter_count(double coeff, unsigned char *paulis, void *op_data)
{
	(void)coeff;
	(void)paulis;

	size_t *count = op_data;
	(*count)++;

	return 0;
}

static int iter_count_onlytwo(
	double coeff, unsigned char *paulis, void *op_data)
{
	(void)coeff;
	(void)paulis;

	size_t *count = op_data;
	(*count)++;
	if (*count == 2)
		return 77;

	return 0;
}

struct iter_store {
	size_t idx;
	size_t num_qubits;
	double coeffs[128];
	unsigned char paulis[128];
};

static int iter_store(double coeff, unsigned char *paulis, void *op_data)
{
	struct iter_store *is = op_data;

	is->coeffs[is->idx] = coeff;
	for (size_t i = 0; i < is->num_qubits; i++) {
		size_t store_idx = is->idx * is->num_qubits + i;
		if (store_idx < 128)
			is->paulis[store_idx] = paulis[i];
	}
	is->idx++;

	return 0;
}

static int test_iter0(void)
{
	const struct test_data td = TEST_DATA[0];
	data_id fid = data_open(td.filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("open file: %s", td.filename);
		return -1;
	}

	size_t num_qubits, num_terms;
	if (data_hamil_getnums(fid, &num_qubits, &num_terms) < 0) {
		TEST_FAIL("wron number of qubits and terms: %zu, %zu",
			num_qubits, num_terms);
		goto err;
	}

	size_t count = 0;
	if (data_hamil_foreach(fid, iter_count, &count) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	if (count != td.num_terms) {
		TEST_FAIL("number of iterations: %zu", count);
		goto err;
	}

	count = 0;
	if (data_hamil_foreach(fid, iter_count_onlytwo, &count) == 0) {
		TEST_FAIL("iteration didn't terminate early");
		goto err;
	}
	if (count != 2) {
		TEST_FAIL("number of iterations: %zu", count);
		goto err;
	}
	data_close(fid);
	return 0;
err:
	data_close(fid);
	return -1;
}

static int test_iter1(void)
{
	const struct test_data td = TEST_DATA[1];
	data_id fid = data_open(td.filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("open file: %s", td.filename);
		return -1;
	}

	size_t num_qubits, num_terms;
	if (data_hamil_getnums(fid, &num_qubits, &num_terms) < 0) {
		TEST_FAIL("wron number of qubits and terms: %zu, %zu",
			num_qubits, num_terms);
		goto err;
	}

	size_t count = 0;
	if (data_hamil_foreach(fid, iter_count, &count) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	if (count != td.num_terms) {
		TEST_FAIL("number of iterations: %zu", count);
		goto err;
	}

	count = 0;
	if (data_hamil_foreach(fid, iter_count_onlytwo, &count) == 0) {
		TEST_FAIL("iteration didn't terminate early");
		goto err;
	}
	if (count != 2) {
		TEST_FAIL("number of iterations: %zu", count);
		goto err;
	}

	struct iter_store is;
	is.idx = 0;
	is.num_qubits = num_qubits;
	if (data_hamil_foreach(fid, iter_store, &is) != 0) {
		TEST_FAIL("iteration terminated early");
		goto err;
	}
	double exp_coeff[] = { 0.64604963, 0.16592673, 0.90327525, -0.18683327,
		-0.18315831, 0.57830137, 0.71210119, -0.96550733, 0.21017606,
		-0.84378561 };
	unsigned char exp_paulis[] = { 0, 3, 1, 3, 0, 3, 0, 1, 0, 1, 1, 3, 0, 3,
		0, 2, 1, 3, 1, 0, 0, 3, 2, 1, 2, 3, 1, 1, 1, 1 };
	for (size_t i = 0; i < num_terms; i++) {
		if (fabs(is.coeffs[i] - exp_coeff[i]) > MARGIN) {
			TEST_FAIL("coeff[%zu] stored: %f vs. expected: %f", i,
				is.coeffs[i], exp_coeff[i]);
			goto err;
		}
	}
	for (size_t i = 0; i < num_terms * num_qubits; i++) {
		if (is.paulis[i] != exp_paulis[i]) {
			TEST_FAIL("pauli[%zu] stored: %u vs. expected: %u", i,
				is.paulis[i], exp_paulis[i]);
			goto err;
		}
	}

	data_close(fid);
	return 0;
err:
	data_close(fid);
	return -1;
}

static int test_iter(void)
{
	if (test_iter0() < 0)
		return -1;
	if (test_iter1() < 0)
		return -1;

	return 0;
}

static void TEST_MAIN(void)
{
	if (test_getnums() < 0) {
		TEST_FAIL("getnums");
		return;
	}
	if (test_getnorm() < 0) {
		TEST_FAIL("getnorm");
		return;
	}
	if (test_iter() < 0) {
		TEST_FAIL("iter");
		return;
	}
}
