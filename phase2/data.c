#include "c23_compat.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "hdf5.h"

#include "phase2/data.h"
#include "phase2/world.h"
#include <complex.h>

static struct world WD;

data_id data_open(const char *filename)
{
	if (world_info(&WD) != WORLD_READY)
		return DATA_INVALID_FID;
	const hid_t acc_plist = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(acc_plist, MPI_COMM_WORLD, MPI_INFO_NULL);

	const hid_t fid = H5Fopen(filename, H5F_ACC_RDWR, acc_plist);
	if (fid == H5I_INVALID_HID)
		return DATA_INVALID_FID;

	return fid;
}

void data_close(const data_id fid)
{
	H5Fclose(fid);
}

int data_grp_create(data_id fid, const char *grp_name)
{
	int rt = -1;

	const hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
	if (lcpl == H5I_INVALID_HID)
		goto ex_lcpl;
	// ensure that intermediate groups are created automatically
	if (H5Pset_create_intermediate_group(lcpl, 1) < 0)
		goto ex_prop;
	// use UTF-8 encoding for link names
	if (H5Pset_char_encoding(lcpl, H5T_CSET_UTF8) < 0)
		goto ex_prop;
	const hid_t group =
		H5Gcreate(fid, grp_name, lcpl, H5P_DEFAULT, H5P_DEFAULT);
	if (group == H5I_INVALID_HID)
		goto ex_group;

	rt = 0;
	H5Gclose(group);
ex_group:
ex_prop:
	H5Pclose(lcpl);
ex_lcpl:
	return rt;
}

#define DEFINE_DATA_ATTR_READ(suff, type, h5_type)                             \
	int data_attr_read_##suff(data_id fid, const char *grp_name,           \
		const char *attr_name, type *a)                                \
	{                                                                      \
		int rt = -1;                                                   \
                                                                               \
		const hid_t grp_id = H5Gopen(fid, grp_name, H5P_DEFAULT);      \
		if (grp_id == H5I_INVALID_HID)                                 \
			goto ex_grp_id;                                        \
		const hid_t attr_id = H5Aopen(grp_id, attr_name, H5P_DEFAULT); \
		if (attr_id == H5I_INVALID_HID)                                \
			goto ex_attr_id;                                       \
		if (H5Aread(attr_id, h5_type, a) < 0)                          \
			goto ex_h5aread;                                       \
                                                                               \
		rt = 0;                                                        \
ex_h5aread:                                                                    \
		H5Aclose(attr_id);                                             \
ex_attr_id:                                                                    \
		H5Gclose(grp_id);                                              \
ex_grp_id:                                                                     \
		return rt;                                                     \
	}

DEFINE_DATA_ATTR_READ(i, int, H5T_NATIVE_INT);
DEFINE_DATA_ATTR_READ(ul, unsigned long, H5T_NATIVE_ULONG);
DEFINE_DATA_ATTR_READ(dbl, double, H5T_NATIVE_DOUBLE);

#define DEFINE_DATA_ATTR_WRITE(suff, type, h5_type)                            \
	int data_attr_write_##suff(data_id fid, const char *grp_name,          \
		const char *attr_name, type a)                                 \
	{                                                                      \
		int rt = -1;                                                   \
                                                                               \
		const hid_t grp_id = H5Gopen(fid, grp_name, H5P_DEFAULT);      \
		if (grp_id == H5I_INVALID_HID)                                 \
			goto ex_group;                                         \
		const hid_t acpl = H5Pcreate(H5P_ATTRIBUTE_CREATE);            \
		if (acpl == H5I_INVALID_HID)                                   \
			goto ex_acpl;                                          \
		if (H5Pset_char_encoding(acpl, H5T_CSET_UTF8) < 0)             \
			goto ex_fspace;                                        \
		const hid_t fspace = H5Screate(H5S_SCALAR);                    \
		if (fspace == H5I_INVALID_HID)                                 \
			goto ex_fspace;                                        \
		const hid_t attr_id = H5Acreate2(grp_id, attr_name, h5_type,   \
			fspace, acpl, H5P_DEFAULT);                            \
		if (attr_id == H5I_INVALID_HID)                                \
			goto ex_attr;                                          \
		if (H5Awrite(attr_id, h5_type, &a) < 0)                        \
			goto ex_write;                                         \
                                                                               \
		rt = 0;                                                        \
ex_write:                                                                      \
		H5Aclose(attr_id);                                             \
ex_attr:                                                                       \
		H5Sclose(fspace);                                              \
ex_fspace:                                                                     \
		H5Pclose(acpl);                                                \
ex_acpl:                                                                       \
		H5Gclose(grp_id);                                              \
ex_group:                                                                      \
		return rt;                                                     \
	}

DEFINE_DATA_ATTR_WRITE(i, int, H5T_NATIVE_INT);
DEFINE_DATA_ATTR_WRITE(ul, unsigned long, H5T_NATIVE_ULONG);
DEFINE_DATA_ATTR_WRITE(dbl, double, H5T_NATIVE_DOUBLE);

struct muldet_handle {
	hid_t stprep_grp_id;
	hid_t muldet_grp_id;
};

static int multidet_open(const data_id fid, struct muldet_handle *md)
{
	const hid_t sp_id = H5Gopen(fid, DATA_STPREP, H5P_DEFAULT);
	if (sp_id == H5I_INVALID_HID)
		goto err_stprep;
	const hid_t md_id = H5Gopen(sp_id, DATA_STPREP_MULTIDET, H5P_DEFAULT);
	if (md_id == H5I_INVALID_HID)
		goto err_muldet;

	md->stprep_grp_id = sp_id;
	md->muldet_grp_id = md_id;

	return 0;

err_muldet:
	H5Gclose(sp_id);
err_stprep:
	return -1;
}

static void multidet_close(struct muldet_handle md)
{
	H5Gclose(md.muldet_grp_id);
	H5Gclose(md.stprep_grp_id);
}

int data_multidet_getnums(data_id fid, uint32_t *nqb, size_t *ndets)
{
	int rt = -1;

	struct muldet_handle md;
	if (multidet_open(fid, &md) < 0)
		return -1;

	const hid_t dset_id = H5Dopen2(
		md.muldet_grp_id, DATA_STPREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID)
		goto ex_md_open;

	const hid_t dsp_id = H5Dget_space(dset_id);
	hsize_t dsp_dims[2];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 2)
		goto ex_dims;

	*ndets = dsp_dims[0];
	*nqb = (uint32_t)dsp_dims[1];

	rt = 0;
ex_dims:
	H5Sclose(dsp_id);
	H5Dclose(dset_id);
ex_md_open:
	multidet_close(md);
	return rt;
}

static int multidet_read_data(data_id fid, double *cfs, unsigned char *dets)
{
	int rt = -1;

	struct muldet_handle md;
	if (multidet_open(fid, &md) < 0)
		goto ex_multidet_open;

	const hid_t dset_cfs_id = H5Dopen2(
		md.muldet_grp_id, DATA_STPREP_MULTIDET_COEFFS, H5P_DEFAULT);
	if (dset_cfs_id == H5I_INVALID_HID)
		goto ex_coeffs_open;
	if (H5Dread(dset_cfs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, cfs) < 0)
		goto ex_coeffs_read;

	const hid_t dset_dets_id = H5Dopen2(
		md.muldet_grp_id, DATA_STPREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_dets_id == H5I_INVALID_HID)
		goto ex_dets_open;
	if (H5Dread(dset_dets_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, dets) < 0)
		goto ex_dets_read;

	rt = 0;

ex_dets_read:
	H5Dclose(dset_dets_id);
ex_dets_open:
ex_coeffs_read:
	H5Dclose(dset_cfs_id);
ex_coeffs_open:
	multidet_close(md);
ex_multidet_open:
	return rt;
}

int data_multidet_foreach(data_id fid,
	int (*op)(_Complex double cf, uint64_t idx, void *), void *op_data)
{
	int rt = -1, rc = 0;
	uint32_t nqb;
	size_t ndets;

	if (data_multidet_getnums(fid, &nqb, &ndets) < 0)
		goto ex_getnums;

	double *cfs = malloc(sizeof *cfs * 2 * ndets);
	if (!cfs)
		goto ex_alloc_coeffs;
	unsigned char *dets = malloc(sizeof *dets * ndets * nqb);
	if (!dets)
		goto ex_alloc_dets;

	/* Read the content of the data file */
	if (multidet_read_data(fid, cfs, dets) < 0)
		goto ex_data_read;

	for (size_t i = 0; i < ndets; i++) {
		uint64_t idx = 0;
		for (size_t j = 0; j < nqb; j++) {
			idx += dets[i * nqb + j] << j;
		}
		_Complex double cf = CMPLX(cfs[2 * i], cfs[2 * i + 1]);
		rc = op(cf, idx, op_data);
		/* This isn't an error, but rather the user telling us to
		   short-circuit the iteration. */
		if (rc != 0)
			break;
	}
	rt = rc;

ex_data_read:
	free(dets);
ex_alloc_dets:
	free(cfs);
ex_alloc_coeffs:
ex_getnums:
	return rt;
}

int data_hamil_getnums(data_id fid, uint32_t *nqb, size_t *nterms)
{
	int rt = -1;

	const hid_t grp_id = H5Gopen(fid, DATA_HAMIL, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID)
		goto ex_grp;
	const hid_t dset_id = H5Dopen2(grp_id, DATA_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID)
		goto ex_dset;
	const hid_t dsp_id = H5Dget_space(dset_id);
	hsize_t dsp_dims[2];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 2)
		goto ex_dims;

	*nterms = dsp_dims[0];
	*nqb = (uint32_t)dsp_dims[1];

	rt = 0;
ex_dims:
	H5Sclose(dsp_id);
	H5Dclose(dset_id);
ex_dset:
	H5Gclose(grp_id);
ex_grp:
	return rt;
}

int data_hamil_getnorm(data_id fid, double *norm)
{
	return data_attr_read_dbl(fid, DATA_HAMIL, DATA_HAMIL_NORM, norm);
}

static int hamil_read_data(data_id fid, double *cfs, unsigned char *paulis)
{
	int rt = -1;

	const hid_t grp_id = H5Gopen(fid, DATA_HAMIL, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID)
		goto ex_grp;
	const hid_t dset_coeffs_id =
		H5Dopen2(grp_id, DATA_HAMIL_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID)
		goto ex_coeffs;
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, cfs) < 0)
		goto ex_coeffs_read;

	const hid_t dset_paulis_id =
		H5Dopen2(grp_id, DATA_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_paulis_id == H5I_INVALID_HID)
		goto ex_paulis;
	if (H5Dread(dset_paulis_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, paulis) < 0)
		goto ex_pauli_read;

	rt = 0;
ex_pauli_read:
	H5Dclose(dset_paulis_id);
ex_paulis:
ex_coeffs_read:
	H5Dclose(dset_coeffs_id);
ex_coeffs:
	H5Gclose(grp_id);
ex_grp:
	return rt;
}

int data_hamil_foreach(const data_id fid,
	int (*op)(double, unsigned char *, void *), void *op_data)
{
	int rt = -1, rc = 0;
	uint32_t nqb;
	size_t nterms;

	if (data_hamil_getnums(fid, &nqb, &nterms) < 0)
		return -1;

	double *cfs = malloc(sizeof *cfs * nterms);
	if (!cfs)
		goto ex_coeffs_alloc;

	unsigned char *paulis = malloc(sizeof *paulis * nqb * nterms);
	if (!paulis)
		goto ex_paulis_alloc;
	if (hamil_read_data(fid, cfs, paulis) < 0)
		goto ex_hamil_read;
	unsigned char *paustr = malloc(sizeof *paustr * nqb);
	if (!paustr)
		goto ex_paustr_alloc;

	for (size_t i = 0; i < nterms; i++) {
		for (size_t j = 0; j < nqb; j++)
			paustr[j] = paulis[i * nqb + j];
		rc = op(cfs[i], paustr, op_data);
		if (rc != 0)
			break;
	}

	rt = rc;
	free(paustr);
ex_paustr_alloc:
ex_hamil_read:
	free(paulis);
ex_paulis_alloc:
	free(cfs);
ex_coeffs_alloc:
	return rt;
}

int data_res_write(data_id fid, const char *grp_name, const char *dset_name,
	const _Complex double *vals, const size_t nvals)
{
	int rt = -1;

	const hid_t grp_id = H5Gopen(fid, grp_name, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID)
		goto ex_open;
	const hid_t dspace_id =
		H5Screate_simple(2, (hsize_t[]){ nvals, 2 }, NULL);
	if (dspace_id == H5I_INVALID_HID)
		goto ex_fspace;

	const hid_t dset_id = H5Dcreate2(grp_id, dset_name, H5T_IEEE_F64LE,
		dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID)
		goto ex_dset;
	if (H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, dspace_id,
		    H5P_DEFAULT, vals) < 0)
		goto ex_dset_write;

	rt = 0;
ex_dset_write:
	H5Dclose(dset_id);
ex_dset:
	H5Sclose(dspace_id);
ex_fspace:
	H5Gclose(grp_id);
ex_open:
	return rt;
}

int data_circ_trott_read_values_test(
	const data_id fid, double *vals[2], const size_t nvals)
{
	int rt = -1;

	double *v = malloc(sizeof(double) * 2 * nvals);
	if (!v)
		goto ex_malloc;

	const hid_t grp_id = H5Gopen(fid, DATA_CIRCTROTT, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID)
		goto ex_open;
	const hid_t dset = H5Dopen2(grp_id, DATA_CIRCTROTT_VALUES, H5P_DEFAULT);
	if (dset == H5I_INVALID_HID)
		goto ex_dset;
	/* clang-format off */
	if (H5Dread(dset, H5T_NATIVE_DOUBLE,
		H5S_ALL, H5S_ALL, H5P_DEFAULT, v) < 0)
		goto ex_read;
	/* clang-format on */

	for (size_t i = 0; i < nvals; i++) {
		vals[0][i] = v[2 * i];
		vals[1][i] = v[2 * i + 1];
	}

	rt = 0;
ex_read:
	H5Dclose(dset);
ex_dset:
	H5Gclose(grp_id);
ex_open:
	free(v);
ex_malloc:
	return rt;
}