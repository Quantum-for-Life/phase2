#include <complex.h>
#include <stdio.h>

#include "data.h"

#include "test.h"
#include "test_data.h"

#define MARGIN (10e-7)

static int
test_get_nums(void)
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

		size_t num_qubits = 0, num_terms = 0;
		if (data2_hamil_getnums(fid, &num_qubits, &num_terms) < 0) {
			TEST_FAIL("read multidet getnums()");
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
		data2_close(fid);
	}

	return rc;
}

int
test_data_hamil()
{
	if (test_get_nums() < 0)
		goto err;

exit:
	return 0;
err:
	return -1;
}
