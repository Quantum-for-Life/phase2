#include "data.h"

#include "test.h"

#ifndef TEST_DATA_SRC
#error
#endif

#define NUM_TEST_FILES (2)
const struct test_data {
	const char *filename;
	size_t	    num_qubits;
	size_t	    num_terms;
	size_t	    num_dets;
} TEST_DATA[NUM_TEST_FILES] = {
	{ .filename	    = TEST_DATA_SRC "/H2O_CAS56.h5",
		.num_qubits = 10,
		.num_terms  = 251,
		.num_dets   = 1 },
	{ .filename	    = TEST_DATA_SRC "/case-d9f603dc.h5_solved",
		.num_qubits = 3,
		.num_terms  = 10,
		.num_dets   = 3 }
};

int
test_data_open(const char *filename);

int
test_data_open_nonexist(void);

int test_data_multidet(data_id, struct test_data);

int
main(void)
{
	/* Check open/close facility */
	if (test_data_open(TEST_DATA[0].filename) < 0)
		goto error;
	if (test_data_open_nonexist() < 0)
		goto error;

	/* Check if num_qubits, num_dets is retrieved correctly */
	for (size_t i = 0; i < NUM_TEST_FILES; i++) {
		data_id	    fid;
		const char *filename = TEST_DATA[i].filename;
		if ((fid = data2_open(filename)) == DATA_INVALID_FID) {
			TEST_FAIL("open test file: %s", filename);
			goto error;
		}
		if (test_data_multidet(fid, TEST_DATA[i]) < 0) {
			TEST_FAIL("multidet test file: %s", filename);
			goto error;
		}
		data2_close(fid);
	}

	/* test iteration */

	return 0;
error:
	return -1;
}