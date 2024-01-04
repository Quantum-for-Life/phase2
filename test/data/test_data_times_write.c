/*
 * Test routine: data2_times_update() from data2.h by creating a fresh H5 file
 * and updating a time series specified in that file.  Check if the values are
 * written correctly.  Remove temporary file.
 */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include <hdf5.h>

#include <test.h>

#define SIZE (5)
#define MARGIN (10e-6)

static struct {
	double		t;
	_Complex double v;
} tst_ts[SIZE] = {
	{ 1.1,  0.0 + _Complex_I * 000.01 },
	{ 13.1, 1.8 + _Complex_I * 000.11 },
	{ 0.42, 2.7 + _Complex_I * 001.11 },
	{ 24,   3.6 + _Complex_I * 011.11 },
	{ 0.01, 4.5 + _Complex_I * 111.11 },
};

static char filename[L_tmpnam];

int
prepare_test_file(hid_t file_id)
{
	return 0;
}

int
test_data_times_write(void)
{
	int   rc = 0;
	hid_t file_id;

	/* Create a temporary H5 file */
	if (!tmpnam(filename)) {
		TEST_FAIL("generate temp filename");
		return -1;
	}

	file_id = H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		TEST_FAIL("create H5 file");
		goto error;
	}

	if (prepare_test_file(file_id) < 0) {
		TEST_FAIL("prepare test file");
		goto error;
	}

	rc = 0;
	goto exit;
error:
	rc = -1;
exit:
	/* Delete temporary file */
	H5Fclose(file_id);
	if (remove(filename) != 0) {
		TEST_FAIL("remove temp file");
		return -1;
	}

	return rc;
}