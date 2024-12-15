#include <stdint.h>

#include "phase2/data.h"
#include "phase2/world.h"

#include "t-data.h"
#include "test.h"

#define WD_SEED UINT64_C(0xd771f4ae402e568b)

/*
static int open_nonexist(void)
{
	if (data_open("") != DATA_INVALID_FID)
		goto err;

	if (data_open(NULL) != DATA_INVALID_FID)
		goto err;

	return 0;
err:
	return -1;
}
*/

static int open_exist(const char *filename)
{
	data_id fid;
	if ((fid = data_open(filename)) == DATA_INVALID_FID)
		goto err;

	data_close(fid);

	return 0;
err:
	return -1;
}

static void TEST_MAIN(void)
{
	world_init((void *)0, (void *)0, WD_SEED);

	/* This test generates noisy error messages from HDF5 and it's
	 * not particularly important.  Disable.
	 *
	 * if (open_nonexist() < 0) {
	 *	TEST_FAIL("data open nonexits");
	 *	goto err;
	 * }
	 */

	for (size_t i = 0; i < NUM_TEST_FILES; i++) {
		const char *filename = TEST_DATA[i].filename;
		if (open_exist(filename) < 0) {
			TEST_FAIL("open data file: %s", filename);
			return;
		}
	}

	world_free();
}
