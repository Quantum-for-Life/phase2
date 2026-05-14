/*
 * Round-trip the /circ_trott group via the per-step write
 * API: data_circ_init pre-allocates a NaN-padded values
 * dataset; data_circ_write_step hyperslab-writes one row
 * at a time.  Reopen, read back via data_attr_read (delta)
 * and direct H5Dread on rank 0 (the values dataset);
 * confirm bit-for-bit match.
 */
#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include <hdf5.h>

#include "phase2/data.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x3a1e27cce387af66)

#define SIZE (5)
#define MARGIN (1.0e-6)

_Complex double tst_vals[SIZE] = {
	CMPLX(0.0, 000.01),
	CMPLX(1.8, 000.11),
	CMPLX(2.7, 001.11),
	CMPLX(3.6, 011.11),
	CMPLX(4.5, 111.11),
};

static double delta = 0.3224;

static char *FILENAME = "/tmp/G1w1Clar2ZLovBir2cGYUbCxgIaV4";

int main(void)
{
	struct world_info wd;
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&wd);

	/* Rank 0 creates an empty file; data_circ_init creates
	 * the /circ_trott group + values dataset. */
	if (wd.rank == 0) {
		remove(FILENAME);
		const hid_t file_id = H5Fcreate(FILENAME, H5F_ACC_EXCL,
			H5P_DEFAULT, H5P_DEFAULT);
		if (file_id == H5I_INVALID_HID)
			TEST_FAIL("create H5 file");
		H5Fclose(file_id);
	}

	data_id fid = data_open(FILENAME);
	if (fid == DATA_INVALID_FID) {
		TEST_FAIL("data: open file");
		goto ex_dat2_open;
	}

	if (data_circ_init(fid, DATA_CIRCTROTT, SIZE) < 0)
		TEST_FAIL("data_circ_init");
	if (data_attr_write(
		    fid, DATA_CIRCTROTT, DATA_CIRCTROTT_DELTA, delta) < 0)
		TEST_FAIL("data_attr_write delta");
	for (size_t i = 0; i < SIZE; i++)
		if (data_circ_write_step(
			    fid, DATA_CIRCTROTT, i, tst_vals[i]) < 0)
			TEST_FAIL("data_circ_write_step %zu", i);
	data_close(fid);

	/* Read delta back through the collective data_attr_read
	 * path; rank 0 reads, all ranks Bcast. */
	fid = data_open(FILENAME);
	if (fid == DATA_INVALID_FID)
		TEST_FAIL("reopen for delta read");
	double d;
	if (data_attr_read(fid, DATA_CIRCTROTT, DATA_CIRCTROTT_DELTA, &d) < 0)
		TEST_FAIL("read delta");
	if (fabs(d - delta) > MARGIN)
		TEST_FAIL("wrong value of delta: %f", d);
	data_close(fid);

	/* Read the values dataset directly on rank 0 to confirm
	 * data_circ_write_step persisted every row. */
	if (wd.rank == 0) {
		const hid_t file_id = H5Fopen(
			FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
		if (file_id == H5I_INVALID_HID)
			TEST_FAIL("open H5 file for readback");
		const hid_t grp_id = H5Gopen(
			file_id, DATA_CIRCTROTT, H5P_DEFAULT);
		const hid_t dset = H5Dopen2(
			grp_id, DATA_CIRCTROTT_VALUES, H5P_DEFAULT);
		_Complex double val_read[SIZE];
		if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
			    H5P_DEFAULT, val_read) < 0)
			TEST_FAIL("H5Dread values");
		for (size_t i = 0; i < SIZE; i++) {
			if (cabs(tst_vals[i] - val_read[i]) > MARGIN)
				TEST_FAIL("wrong value: %f+%fi"
					  " (expected: %f+%fi)",
					creal(val_read[i]), cimag(val_read[i]),
					creal(tst_vals[i]),
					cimag(tst_vals[i]));
		}
		H5Dclose(dset);
		H5Gclose(grp_id);
		H5Fclose(file_id);
	}

ex_dat2_open:
	/* Delete temporary file */
	if (wd.rank == 0 && remove(FILENAME) != 0)
		TEST_FAIL("remove temp file");
	world_free();
}
