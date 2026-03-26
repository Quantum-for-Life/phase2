/*
 * Test routine: data attribute read/write roundtrip.
 *
 * Write int, unsigned long, and double attributes to an HDF5 file via
 * data_attr_write(), close, reopen, read back via data_attr_read(),
 * and verify the values match.
 */
#include "c23_compat.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include <hdf5.h>

#include "phase2/data.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0x7b2f4a91d6e803c5)

#define MARGIN (1.0e-15)

static char *FILENAME = "/tmp/N3kF8xVmPq2AttrTest";

static const char *GRP_NAME = "test_grp";
static const char *ATTR_INT = "val_int";
static const char *ATTR_UL = "val_ul";
static const char *ATTR_DBL = "val_dbl";

static const int VAL_INT = 42;
static const unsigned long VAL_UL = 123456789UL;
static const double VAL_DBL = 3.14159265358979;

int main(void)
{
	struct world_info wd;
	world_init(nullptr, nullptr, WD_SEED);
	world_info(&wd);

	if (wd.rank == 0)
		remove(FILENAME);

	/*
	 * Set up file access property list with parallel I/O access
	 */
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

	hid_t file_id =
		H5Fcreate(FILENAME, H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
	if (file_id == H5I_INVALID_HID)
		TEST_FAIL("create H5 file");

	H5Pclose(plist_id);

	/* Create group and write attributes */
	data_id fid = (data_id)file_id;

	if (data_grp_create(fid, GRP_NAME) < 0)
		TEST_FAIL("data_grp_create");

	if (data_attr_write(fid, GRP_NAME, ATTR_INT, VAL_INT) < 0)
		TEST_FAIL("data_attr_write int");

	if (data_attr_write(fid, GRP_NAME, ATTR_UL, VAL_UL) < 0)
		TEST_FAIL("data_attr_write unsigned long");

	if (data_attr_write(fid, GRP_NAME, ATTR_DBL, VAL_DBL) < 0)
		TEST_FAIL("data_attr_write double");

	H5Fclose(file_id);

	/* Reopen with data_open and read back */
	fid = data_open(FILENAME);
	if (fid == DATA_INVALID_FID)
		TEST_FAIL("data_open");

	int rd_int;
	unsigned long rd_ul;
	double rd_dbl;

	if (data_attr_read(fid, GRP_NAME, ATTR_INT, &rd_int) < 0)
		TEST_FAIL("data_attr_read int");
	TEST_ASSERT(rd_int == VAL_INT,
		"int mismatch: got %d, expected %d", rd_int, VAL_INT);

	if (data_attr_read(fid, GRP_NAME, ATTR_UL, &rd_ul) < 0)
		TEST_FAIL("data_attr_read unsigned long");
	TEST_ASSERT(rd_ul == VAL_UL,
		"unsigned long mismatch: got %lu, expected %lu",
		rd_ul, VAL_UL);

	if (data_attr_read(fid, GRP_NAME, ATTR_DBL, &rd_dbl) < 0)
		TEST_FAIL("data_attr_read double");
	TEST_ASSERT(fabs(rd_dbl - VAL_DBL) < MARGIN,
		"double mismatch: got %.17g, expected %.17g",
		rd_dbl, VAL_DBL);

	data_close(fid);

	/* Delete temporary file */
	if (wd.rank == 0 && remove(FILENAME) != 0)
		TEST_FAIL("remove temp file");

	world_free();
}
