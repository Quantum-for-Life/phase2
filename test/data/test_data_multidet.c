#include "data.h"

#include "test.h"

struct test_data {
	const char *filename;
	size_t      num_qubits;
	size_t      num_terms;
	size_t      num_dets;
};

static int
test_nums(data_id fid, struct test_data td)
{
	size_t num_qubits = 0, num_dets = 0;

	if (data2_multidet_getnums(fid, &num_qubits, &num_dets) < 0) {
		TEST_FAIL("read multidet getnums()");
		goto error;
	}
	if (num_qubits != td.num_qubits) {
		TEST_FAIL("wrong number of qubits: %zu", num_qubits);
		goto error;
	}
	if (num_dets != td.num_dets) {
		TEST_FAIL("wrong number of dets: %zu", num_dets);
		goto error;
	}

	return 0;
error:
	return -1;
}

int
test_data_multidet(data_id fid, struct test_data td)
{
	if (test_nums(fid, td) < 0)
		goto error;

	return 0;
error:
	return -1;
}
