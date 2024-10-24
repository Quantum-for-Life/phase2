#include "c23_compat.h"
#include <stddef.h>
#include <stdint.h>

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

static int data_group_open(const hid_t fid, hid_t *grp_id, const char *name)
{
	const hid_t hid = H5Gopen2(fid, name, H5P_DEFAULT);
	if (hid == H5I_INVALID_HID)
		return -1;
	*grp_id = hid;

	return 0;
}

static void data_group_close(const hid_t grp_id)
{
	H5Gclose(grp_id);
}

static int data_read_attr(hid_t fid, const char *name, hid_t type_id, void *buf)
{
	int rt = -1;

	const hid_t attr_id = H5Aopen(fid, name, H5P_DEFAULT);
	if (attr_id == H5I_INVALID_HID)
		goto ex_attr_id;
	if (H5Aread(attr_id, type_id, buf) < 0)
		goto ex_h5aread;

	rt = 0;
ex_h5aread:
	H5Aclose(attr_id);
ex_attr_id:
	return rt;
}

struct muldet_handle {
	hid_t stprep_grp_id;
	hid_t muldet_grp_id;
};

static int multidet_open(const data_id fid, struct muldet_handle *md)
{
	hid_t sp_id, md_id;

	if (data_group_open(fid, &sp_id, DATA_STPREP) < 0)
		goto err_grp_open_stprep;
	if (data_group_open(sp_id, &md_id, DATA_STPREP_MULTIDET) < 0)
		goto err_grp_open_muldet;

	md->stprep_grp_id = sp_id;
	md->muldet_grp_id = md_id;

	return 0;

err_grp_open_muldet:
	data_group_close(sp_id);
err_grp_open_stprep:
	return -1;
}

static void multidet_close(struct muldet_handle md)
{
	data_group_close(md.muldet_grp_id);
	data_group_close(md.stprep_grp_id);
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

	hid_t grp_id;
	if (data_group_open(fid, &grp_id, DATA_HAMIL) < 0)
		goto ex_open;

	const hid_t dset_id = H5Dopen2(grp_id, DATA_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID)
		goto ex_read;

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
ex_read:
	data_group_close(grp_id);
ex_open:
	return rt;
}

int data_hamil_getnorm(data_id fid, double *norm)
{
	int rt = -1;

	hid_t grp_id;
	if (data_group_open(fid, &grp_id, DATA_HAMIL) < 0)
		goto ex_open;

	const hid_t attr_norm_id =
		H5Aopen(grp_id, DATA_HAMIL_NORM, H5P_DEFAULT);
	if (attr_norm_id == H5I_INVALID_HID)
		goto ex_attr_open;

	double n;
	if (H5Aread(attr_norm_id, H5T_NATIVE_DOUBLE, &n) < 0)
		goto ex_attr_read;
	*norm = n;

	rt = 0;
ex_attr_read:
	H5Aclose(attr_norm_id);
ex_attr_open:
	data_group_close(grp_id);
ex_open:
	return rt;
}

static int hamil_read_data(data_id fid, double *cfs, unsigned char *paulis)
{
	int rt = -1;

	hid_t grp_id;
	if (data_group_open(fid, &grp_id, DATA_HAMIL) < 0)
		goto ex_data_group_open;

	const hid_t dset_coeffs_id =
		H5Dopen2(grp_id, DATA_HAMIL_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID)
		goto ex_coeffs_open;

	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, cfs) < 0)
		goto ex_coeffs_read;

	const hid_t dset_paulis_id =
		H5Dopen2(grp_id, DATA_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_paulis_id == H5I_INVALID_HID)
		goto ex_paulis_open;

	if (H5Dread(dset_paulis_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, paulis) < 0)
		goto ex_pauli_read;

	rt = 0;
ex_pauli_read:
	H5Dclose(dset_paulis_id);
ex_paulis_open:
ex_coeffs_read:
	H5Dclose(dset_coeffs_id);
ex_coeffs_open:
	data_group_close(grp_id);
ex_data_group_open:
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

int data_circ_trott_getttrs(data_id fid, double *delta)
{
	int rt = -1;

	hid_t grp_id;
	if (data_group_open(fid, &grp_id, DATA_CIRCTROTT) < 0)
		goto ex_data_group_open;

	if (data_read_attr(grp_id, DATA_CIRCTROTT_TIMEFACTOR, H5T_NATIVE_DOUBLE,
		    delta) < 0)
		goto ex_data_read_attr;

	rt = 0;
ex_data_read_attr:
	data_group_close(grp_id);
ex_data_group_open:
	return rt;
}

int data_write_vals(data_id fid, char *grp_name, char *dset_name,
	_Complex double *vals, size_t nvals)
{
	int rt = -1;

	hid_t grp_id, dspace_id, dset_id;
	if (data_group_open(fid, &grp_id, grp_name) < 0)
		goto ex_open;
	if ((dspace_id = H5Screate_simple(2, (hsize_t[]){ nvals, 2 }, NULL)) ==
		H5I_INVALID_HID)
		goto ex_fspace;

	if ((dset_id = H5Dcreate2(grp_id, dset_name, H5T_IEEE_F64LE, dspace_id,
		     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) == H5I_INVALID_HID)
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
	data_group_close(grp_id);
ex_open:
	return rt;
}

int data_circ_trott_read_values_test(
	const data_id fid, double *vals[2], const size_t nvals)
{
	int rt = -1;

	hid_t grp_id;
	double *v = malloc(sizeof(double) * 2 * nvals);
	if (!v)
		goto ex_malloc;

	if (data_group_open(fid, &grp_id, DATA_CIRCTROTT) < 0)
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
	data_group_close(grp_id);
ex_open:
	free(v);
ex_malloc:
	return rt;
}

int data_circ_qdrift_getattrs(
	const data_id fid, size_t *nsamples, double *step_size, size_t *depth)
{
	int rt = 0;

	hid_t grp_id;
	if (data_group_open(fid, &grp_id, DATA_CIRCQDRIFT) < 0)
		goto ex_data_group_open;

	if (data_read_attr(grp_id, DATA_CIRCQDRIFT_NUMSAMPLES,
		    H5T_NATIVE_UINT64, nsamples) < 0)
		goto ex_data_read_attr;
	if (data_read_attr(grp_id, DATA_CIRCQDRIFT_STEPSIZE, H5T_NATIVE_DOUBLE,
		    step_size) < 0)
		goto ex_data_read_attr;
	if (data_read_attr(grp_id, DATA_CIRCQDRIFT_DEPTH, H5T_NATIVE_UINT64,
		    depth) < 0)
		goto ex_data_read_attr;

	rt = 0;
ex_data_read_attr:
	data_group_close(grp_id);
ex_data_group_open:
	return rt;
}