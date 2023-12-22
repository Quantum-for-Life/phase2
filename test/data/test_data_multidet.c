#include "data.h"

#include "test.h"

struct test_data {
	const char *filename;
	size_t	    num_qubits;
	size_t	    num_terms;
	size_t	    num_dets;
};

static int
test_num_qubits(data_id fid, struct test_data td)
{
	size_t n = 0;

	if (data2_multidet_num_qubits(fid, &n) < 0) {
		TEST_FAIL("read multidet_num_qubits");
		goto error;
	}
	if (n != td.num_qubits) {
		TEST_FAIL("wrong no of qubits: %zu", n);
		goto error;
	}

	return 0;
error:
	return -1;
}

static int
test_num_terms(data_id fid, struct test_data td)
{
	size_t n = 0;

	if (data2_multidet_num_terms(fid, &n) < 0) {
		TEST_FAIL("read multidet_num_terms");
		goto error;
	}
	if (n != td.num_dets) {
		TEST_FAIL("wrong no. of terms: %zu", n);
		goto error;
	}

	return 0;
error:
	return -1;
}

int
test_data_multidet(data_id fid, struct test_data td)
{
	if (test_num_qubits(fid, td) < 0)
		goto error;
	if (test_num_terms(fid, td) < 0)
		goto error;

	return 0;
error:
	return -1;
}