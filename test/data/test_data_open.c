#include "data.h"

#include "test.h"
#include "test_data.h"

static int
open_nonexist(void)
{
	if (data2_open("") != DATA_INVALID_FID)
		goto error;

	if (data2_open(NULL) != DATA_INVALID_FID)
		goto error;

	return 0;
error:
	return -1;
}

static int
open_exist(const char *filename)
{
	data_id fid;
	if ((fid = data2_open(filename) == DATA_INVALID_FID))
		goto error;

	data2_close(fid);

	return 0;
error:
	return -1;
}

int
test_data_open(void)
{
	data_id fid;

	if (open_nonexist() < 0) {
		TEST_FAIL("data open nonexits");
		goto error;
	}

	for (size_t i = 0; i < NUM_TEST_FILES; i++) {
		const char *filename = TEST_DATA[i].filename;
		if (open_exist(filename) < 0) {
			TEST_FAIL("open data file: %s", filename);
			goto error;
		}
	}

	return 0;
error:
	return -1;
}
