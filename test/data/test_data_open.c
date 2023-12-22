#include "data.h"

int
test_data_open_nonexist(void)
{
	if (data2_open("") != DATA_INVALID_FID)
		goto error;

	if (data2_open(NULL) != DATA_INVALID_FID)
		goto error;

	return 0;
error:
	return -1;
}

int
test_data_open(void)
{
	data_id fid;

	if ((fid = data2_open(TEST_DATA_FILE) == DATA_INVALID_FID))
		goto error;

	data2_close(fid);
	return 0;
error:
	return -1;
}
