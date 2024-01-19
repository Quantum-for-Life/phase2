/*
 * Test routine: data2_times_update() from data2.h by creating a fresh H5 file
 * and updating a time series specified in that file.  Check if the values are
 * written correctly.  Remove temporary file.
 */

#include <complex.h>
#include <data2.h>
#include <stdio.h>
#include <stdlib.h>

#include <hdf5.h>

#include <test.h>

#define SIZE (5)
#define MARGIN (10e-6)

#define H5_GRP_NAME "time_series"
#define H5_GRP_TIMES "times"
#define H5_GRP_VALUES "values"

static struct {
	double		t;
	_Complex double v;
} tst_ts[SIZE] = {
	{1.1,   0.0 + _Complex_I * 000.01},
	{ 13.1, 1.8 + _Complex_I * 000.11},
	{ 0.42, 2.7 + _Complex_I * 001.11},
	{ 24,   3.6 + _Complex_I * 011.11},
	{ 0.01, 4.5 + _Complex_I * 111.11},
};

int iter_update(double *t, _Complex double *v, void *iter_dat)
{
	size_t *i = iter_dat;
	*t	  = tst_ts[*i].t;
	*v	  = tst_ts[*i].v;
	(*i)++;

	return 0;
}

int iter_check(double t, _Complex double v, void *iter_dat)
{
	size_t *i = iter_dat;
	if (t != tst_ts[*i].t || v != tst_ts[*i].v)
		return -1;
	(*i)++;

	return 0;
}

static char filename[L_tmpnam];

enum ret_code {
	ERR = INT_MIN,
	ERR_PREP,
	ERR_CREATEFILE,
	ERR_DELFILE,
	ERR_H5DSET,
	ERR_H5WRITE,
	ERR_DAT2,
	OK = 0,
};

int prepare_dset_times(hid_t grp_id)
{
	enum ret_code rc = OK;

	double times[SIZE];
	for (size_t i = 0; i < SIZE; i++) {
		times[i] = 0.0;
	}

	const hid_t dspace = H5Screate_simple(1, (hsize_t[]){ SIZE }, NULL);
	if (dspace == H5I_INVALID_HID) {
		rc = ERR_H5DSET;
		goto ex_dspace;
	}
	const hid_t dset = H5Dcreate2(grp_id, H5_GRP_TIMES, H5T_NATIVE_DOUBLE,
		dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset == H5I_INVALID_HID) {
		rc = ERR_H5DSET;
		goto ex_dset;
	}

	if (H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, dspace, H5P_DEFAULT,
		    times) < 0) {
		rc = ERR_H5WRITE;
		goto ex_write;
	}

ex_write:
	H5Dclose(dset);
ex_dset:
	H5Sclose(dspace);
ex_dspace:
	return rc;
}

int prepare_dset_values(hid_t grp_id)
{
	enum ret_code rc = OK;

	_Complex double values[SIZE];
	for (size_t i = 0; i < SIZE; i++) {
		values[i] = 0.0;
	}

	const hid_t dspace = H5Screate_simple(2, (hsize_t[]){ SIZE, 2 }, NULL);
	if (dspace == H5I_INVALID_HID) {
		rc = ERR_H5DSET;
		goto ex_dspace;
	}
	const hid_t dset = H5Dcreate2(grp_id, H5_GRP_VALUES, H5T_NATIVE_DOUBLE,
		dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset == H5I_INVALID_HID) {
		rc = ERR_H5DSET;
		goto ex_dset;
	}

	if (H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, dspace, H5P_DEFAULT,
		    values) < 0) {
		rc = ERR_H5WRITE;
		goto ex_write;
	}

ex_write:
	H5Dclose(dset);
ex_dset:
	H5Sclose(dspace);
ex_dspace:
	return rc;
}

int prepare_test_file(hid_t file_id)
{
	enum ret_code rc = OK;

	/* Create main group */
	const hid_t grp_id = H5Gcreate(
		file_id, H5_GRP_NAME, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID) {
		rc = ERR_CREATEFILE;
		goto ex_create;
	}

	/* Write zeroed data sets */
	if (prepare_dset_times(grp_id) < OK) {
		TEST_FAIL("prepare dataset: times");
		rc = ERR_PREP;
		goto ex_prepare;
	}
	if (prepare_dset_values(grp_id) < OK) {
		TEST_FAIL("prepare dataset: values");
		rc = ERR_PREP;
		goto ex_prepare;
	}

ex_prepare:
	H5Gclose(grp_id);
ex_create:
	return rc;
}

int test_data_times_write(void)
{
	enum ret_code rc;

	/* Create a temporary H5 file */
	if (!tmpnam(filename)) {
		TEST_FAIL("generate temp filename");
		rc = ERR_CREATEFILE;
		goto ex_create;
	}
	const hid_t file_id =
		H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		TEST_FAIL("create H5 file");
		rc = ERR_CREATEFILE;
		goto ex_create;
	}
	rc = prepare_test_file(file_id);
	H5Fclose(file_id);
	if (rc < OK) {
		TEST_FAIL("prepare test file");
		rc = ERR_PREP;
		goto ex_prepare;
	}

	const data2_id fid = data2_open(filename);
	if (fid == DATA2_INVALID_FID) {
		TEST_FAIL("data2: reopen file");
		rc = ERR_DAT2;
		goto ex_dat2_open;
	}

	size_t idx = 0;
	if (data2_times_update(fid, iter_update, &idx) < 0) {
		TEST_FAIL("data2: update");
		rc = ERR_DAT2;
		goto ex_dat2_update;
	}

	idx = 0;
	if (data2_times_foreach(fid, iter_check, &idx) < 0) {
		TEST_FAIL("data2: check");
		rc = ERR_DAT2;
	}

ex_dat2_update:
	data2_close(fid);
ex_dat2_open:
ex_prepare:
	/* Delete temporary file */
	// printf("Temp file: %s\n", filename);
	if (remove(filename) != 0) {
		TEST_FAIL("remove temp file");
		rc = ERR_DELFILE;
		goto ex_create;
	}
ex_create:
	return rc == OK ? 0 : -1;
}