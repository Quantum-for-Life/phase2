#include "data.h"

#include "test.h"
#include "test_data.h"

static int
get_nums(data_id fid, struct test_data td)
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

static int
iter_test0(void)
{

	return -1;
}

int
test_data_multidet()
{
	int rc = 0;

	for (size_t i = 0; i < NUM_TEST_FILES; i++) {
		const char *filename = TEST_DATA[i].filename;

		data_id fid = data2_open(filename);
		if (fid == DATA_INVALID_FID) {
			TEST_FAIL("open file: %s", filename);
			rc = -1;
			break;
		}

		if (get_nums(fid, TEST_DATA[i]) < 0) {
			TEST_FAIL("get nums for file: %s", filename);
			rc = -1;
		}

		data2_close(fid);
	}

	if (rc != 0)
		goto error;

	return 0;
error:
	return -1;
}
