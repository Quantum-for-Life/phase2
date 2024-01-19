#include "data2.h"

#include "test.h"
#include "test_data.h"

static int open_nonexist(void)
{
	if (data2_open("") != DATA2_INVALID_FID)
		goto err;

	if (data2_open(NULL) != DATA2_INVALID_FID)
		goto err;

	return 0;
err:
	return -1;
}

static int open_exist(const char *filename)
{
	data2_id fid;
	if ((fid = data2_open(filename) == DATA2_INVALID_FID))
		goto err;

	data2_close(fid);

	return 0;
err:
	return -1;
}

int test_data_open(void)
{
	if (open_nonexist() < 0) {
		TEST_FAIL("data open nonexits");
		goto err;
	}

	for (size_t i = 0; i < NUM_TEST_FILES; i++) {
		const char *filename = TEST_DATA[i].filename;
		if (open_exist(filename) < 0) {
			TEST_FAIL("open data file: %s", filename);
			goto err;
		}
	}

	return 0;
err:
	return -1;
}
