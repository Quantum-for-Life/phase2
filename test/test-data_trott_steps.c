/*
 * Test routine: data_trott_write() from data.h by creating a fresh H5 file
 * Check if the values are written correctly.  Remove temporary file.
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>

#include <hdf5.h>

#include "data.h"
#include "test.h"

#define SIZE (5)
#define MARGIN (1e-6)

#define H5_GRP_NAME "circ_trott"
#define H5_GRP_TIME_FACTOR "time_factor"
#define H5_GRP_VALUES "values"

_Complex double tst_vals[SIZE] = {
	0.0 + _Complex_I * 000.01,
	1.8 + _Complex_I * 000.11,
	2.7 + _Complex_I * 001.11,
	3.6 + _Complex_I * 011.11,
	4.5 + _Complex_I * 111.11,
};

static double time_factor = 0.3224;

static char filename[L_tmpnam];

enum ret_code {
	ERR = INT_MIN,
	ERR_PREP,
	ERR_CREATEFILE,
	ERR_DELFILE,
	ERR_H5ATTR,
	ERR_H5WRITE,
	ERR_H5READ,
	ERR_DAT2,
	ERR_OK = 0
};

static int prepare_test_file(hid_t file_id)
{
	enum ret_code rc = ERR_OK;

	hid_t grp_id, fspace, attr;

	/* Create main group */
	grp_id = H5Gcreate(
		file_id, H5_GRP_NAME, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID) {
		rc = ERR_CREATEFILE;
		goto ex_create;
	}

	// create a scalar (singleton) attribute
	if ((fspace = H5Screate(H5S_SCALAR)) == H5I_INVALID_HID) {
		rc = ERR_H5ATTR;
		goto ex_fspace;
	}
	if ((attr = H5Acreate2(grp_id, H5_GRP_TIME_FACTOR, H5T_IEEE_F64LE,
		     fspace, H5P_DEFAULT, H5P_DEFAULT)) == H5I_INVALID_HID) {
		rc = ERR_H5WRITE;
		goto ex_attr;
	}
	if (H5Awrite(attr, H5T_NATIVE_DOUBLE, &time_factor) < 0) {
		rc = ERR_H5WRITE;
		goto ex_attr_write;
	}

ex_attr_write:
	H5Aclose(attr);
ex_attr:
	H5Sclose(fspace);
ex_fspace:
	H5Gclose(grp_id);
ex_create:
	return rc;
}

int test_data_trott_steps(void)
{
	enum ret_code rc;

	/* Create a temporary H5 file */
	if (!tmpnam(filename)) {
		TEST_FAIL("generate temp filename");
		rc = ERR_CREATEFILE;
		goto ex_create;
	}
	hid_t file_id =
		H5Fcreate(filename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		TEST_FAIL("create H5 file");
		rc = ERR_CREATEFILE;
		goto ex_create;
	}
	rc = prepare_test_file(file_id);
	H5Fclose(file_id);
	if (rc < ERR_OK) {
		TEST_FAIL("prepare test file");
		rc = ERR_PREP;
		goto ex_prepare;
	}

	data_id fid = data_open(filename);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("data: reopen file");
		rc = ERR_DAT2;
		goto ex_dat2_open;
	}

	double tf;
	if (data_circ_trott_getttrs(fid, &tf) < 0) {
		TEST_FAIL("data: read time_factor attribute");
		rc = ERR_DAT2;
		goto ex;
	}
	if (fabs(tf - time_factor) > MARGIN) {
		TEST_FAIL("data: wrong value of the attribute read");
		rc = ERR_DAT2;
		goto ex;
	}

	double tst_vals_re[SIZE], tst_vals_im[SIZE];
	for (size_t i = 0; i < SIZE; i++) {
		tst_vals_re[i] = creal(tst_vals[i]);
		tst_vals_im[i] = cimag(tst_vals[i]);
	}
	data_circ_trott_write_values(
		fid, (double *[]){ tst_vals_re, tst_vals_im }, SIZE);
	data_close(fid);

	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		TEST_FAIL("open H5 file");
		rc = ERR_CREATEFILE;
		goto ex_create;
	}

	hid_t grp_id = H5Gopen(file_id, H5_GRP_NAME, H5P_DEFAULT);
	hid_t dset = H5Dopen2(grp_id, H5_GRP_VALUES, H5P_DEFAULT);
	_Complex double val_read[SIZE];
	if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    val_read) < 0) {
		rc = ERR_H5READ;
		goto ex_readh5;
	}

	for (size_t i = 0; i < SIZE; i++) {
		if (cabs(tst_vals[i] - val_read[i]) > MARGIN) {
			TEST_FAIL("wrong value: %f+%fi (expected: %f+%fi)",
				creal(val_read[i]), cimag(val_read[i]),
				creal(tst_vals[i]), cimag(tst_vals[i]));
			rc = ERR_H5READ;
			goto ex_readh5;
		}
	}

ex_readh5:
	H5Dclose(dset);
	H5Gclose(grp_id);
	H5Fclose(file_id);
	fid = data_open(filename);
ex:
	data_close(fid);
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
	return rc == ERR_OK ? 0 : -1;
}
