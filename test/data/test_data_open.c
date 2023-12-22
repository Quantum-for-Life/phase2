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
test_data_open(const char* filename)
{
	data_id fid;

	if ((fid = data2_open(filename) == DATA_INVALID_FID))
		goto error;

	data2_close(fid);
	return 0;
error:
	return -1;
}
