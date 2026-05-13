#include "c23_compat.h"
#include <errno.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "hdf5.h"

#include "phase2/data.h"
#include "phase2/world.h"
#include <complex.h>

static struct world_info WD;

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
	if (dsp_id == H5I_INVALID_HID)
		goto ex_dsp;
	hsize_t dsp_dims[2];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 2)
		goto ex_dims;

	*ndets = dsp_dims[0];
	*nqb = (uint32_t)dsp_dims[1];

	rt = 0;
ex_dims:
	H5Sclose(dsp_id);
ex_dsp:
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
			idx += (uint64_t)dets[i * nqb + j] << j;
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

/*
 * coeff_matrix readers.  Mirror the multidet pattern (H5Gopen
 * for the subgroup, H5Aread for each attribute, H5Dread with
 * H5T_NATIVE_DOUBLE for the C datasets).
 */

int data_state_prep_kind(const data_id fid, enum stprep_kind *out)
{
	const hid_t sp_id = H5Gopen(fid, DATA_STPREP, H5P_DEFAULT);
	if (sp_id == H5I_INVALID_HID)
		return -ENOENT;

	const htri_t has_md =
		H5Lexists(sp_id, DATA_STPREP_MULTIDET, H5P_DEFAULT);
	const htri_t has_cm =
		H5Lexists(sp_id, DATA_STPREP_COEFFMAT, H5P_DEFAULT);
	H5Gclose(sp_id);

	if (has_md < 0 || has_cm < 0)
		return -EIO;
	if (has_md && has_cm) {
		fprintf(stderr,
			"simul.h5: ambiguous state prep (both "
			"/state_prep/multidet and "
			"/state_prep/coeff_matrix present); "
			"rebuild simul.h5 with exactly one\n");
		return -EINVAL;
	}
	if (!has_md && !has_cm) {
		fprintf(stderr,
			"simul.h5: no state-prep subgroup found "
			"(expected /state_prep/multidet or "
			"/state_prep/coeff_matrix)\n");
		return -ENOENT;
	}

	*out = has_md ? STPREP_MULTIDET : STPREP_COEFF_MATRIX;
	return 0;
}

struct coeffmat_handle {
	hid_t stprep_grp_id;
	hid_t coeffmat_grp_id;
};

static int coeffmat_open(const data_id fid, struct coeffmat_handle *cm)
{
	const hid_t sp_id = H5Gopen(fid, DATA_STPREP, H5P_DEFAULT);
	if (sp_id == H5I_INVALID_HID)
		goto err_stprep;
	const hid_t cm_id = H5Gopen(sp_id, DATA_STPREP_COEFFMAT, H5P_DEFAULT);
	if (cm_id == H5I_INVALID_HID)
		goto err_cm;

	cm->stprep_grp_id = sp_id;
	cm->coeffmat_grp_id = cm_id;
	return 0;

err_cm:
	H5Gclose(sp_id);
err_stprep:
	return -1;
}

static void coeffmat_close(struct coeffmat_handle cm)
{
	H5Gclose(cm.coeffmat_grp_id);
	H5Gclose(cm.stprep_grp_id);
}

static int read_u32_attr(hid_t grp_id, const char *name, uint32_t *out)
{
	const hid_t aid = H5Aopen(grp_id, name, H5P_DEFAULT);
	if (aid == H5I_INVALID_HID)
		return -1;
	const int rt = H5Aread(aid, H5T_NATIVE_UINT32, out) < 0 ? -1 : 0;
	H5Aclose(aid);
	return rt;
}

static int read_u8_attr(hid_t grp_id, const char *name, int *out)
{
	const hid_t aid = H5Aopen(grp_id, name, H5P_DEFAULT);
	if (aid == H5I_INVALID_HID)
		return -1;
	uint8_t v = 0;
	const int rt = H5Aread(aid, H5T_NATIVE_UINT8, &v) < 0 ? -1 : 0;
	H5Aclose(aid);
	if (rt == 0)
		*out = v ? 1 : 0;
	return rt;
}

static int read_double_attr(hid_t grp_id, const char *name, double *out)
{
	const hid_t aid = H5Aopen(grp_id, name, H5P_DEFAULT);
	if (aid == H5I_INVALID_HID)
		return -1;
	const int rt = H5Aread(aid, H5T_NATIVE_DOUBLE, out) < 0 ? -1 : 0;
	H5Aclose(aid);
	return rt;
}

static int validate_C_shape(hid_t grp_id, const char *dset_name,
	uint32_t n_sites, uint32_t n_occ)
{
	const hid_t did = H5Dopen2(grp_id, dset_name, H5P_DEFAULT);
	if (did == H5I_INVALID_HID)
		return -1;
	int rt = -1;
	const hid_t sid = H5Dget_space(did);
	if (sid == H5I_INVALID_HID)
		goto ex_dset;
	hsize_t dims[2] = { 0, 0 };
	if (H5Sget_simple_extent_dims(sid, dims, NULL) != 2)
		goto ex_space;
	if (dims[0] != n_sites || dims[1] != n_occ)
		goto ex_space;
	const hid_t tid = H5Dget_type(did);
	if (tid == H5I_INVALID_HID)
		goto ex_space;
	const H5T_class_t cls = H5Tget_class(tid);
	const size_t sz = H5Tget_size(tid);
	H5Tclose(tid);
	if (cls != H5T_FLOAT || sz != sizeof(double))
		goto ex_space;
	rt = 0;
ex_space:
	H5Sclose(sid);
ex_dset:
	H5Dclose(did);
	return rt;
}

int data_coeff_matrix_getnums(const data_id fid, uint32_t *nqb,
	uint32_t *n_sites, uint32_t *n_alpha, uint32_t *n_beta,
	int *closed_shell, int *tapered)
{
	int rt = -1;
	struct coeffmat_handle cm;
	if (coeffmat_open(fid, &cm) < 0)
		return -1;

	uint32_t v_nqb = 0, v_ns = 0, v_na = 0, v_nb = 0;
	int v_cs = 0, v_tap = 0;

	if (read_u32_attr(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_NQB, &v_nqb)
		< 0)
		goto ex;
	if (read_u32_attr(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_NS, &v_ns)
		< 0)
		goto ex;
	if (read_u32_attr(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_NA, &v_na)
		< 0)
		goto ex;
	if (read_u32_attr(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_NB, &v_nb)
		< 0)
		goto ex;
	if (read_u8_attr(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_CS, &v_cs)
		< 0)
		goto ex;
	if (read_u8_attr(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_TAP, &v_tap)
		< 0)
		goto ex;

	*nqb = v_nqb;
	*n_sites = v_ns;
	*n_alpha = v_na;
	*n_beta = v_nb;
	*closed_shell = v_cs;
	*tapered = v_tap;
	rt = 0;
ex:
	coeffmat_close(cm);
	return rt;
}

static int read_C_dset(hid_t grp_id, const char *name, uint32_t n_sites,
	uint32_t n_occ, double *buf)
{
	if (validate_C_shape(grp_id, name, n_sites, n_occ) < 0)
		return -1;

	const hid_t did = H5Dopen2(grp_id, name, H5P_DEFAULT);
	if (did == H5I_INVALID_HID)
		return -1;
	int rt = -1;
	if (H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf)
		< 0)
		goto ex;
	rt = 0;
ex:
	H5Dclose(did);
	return rt;
}

int data_coeff_matrix_read(
	const data_id fid, double *C_alpha, double *C_beta)
{
	int rt = -1;

	uint32_t nqb, n_sites, n_alpha, n_beta;
	int closed_shell, tapered;
	if (data_coeff_matrix_getnums(fid, &nqb, &n_sites, &n_alpha, &n_beta,
		    &closed_shell, &tapered) < 0)
		return -1;

	if (closed_shell && C_beta != NULL)
		return -1;
	if (!closed_shell && C_beta == NULL)
		return -1;

	struct coeffmat_handle cm;
	if (coeffmat_open(fid, &cm) < 0)
		return -1;

	if (read_C_dset(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_CA, n_sites,
		    n_alpha, C_alpha) < 0)
		goto ex;
	if (!closed_shell) {
		if (read_C_dset(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_CB,
			    n_sites, n_beta, C_beta) < 0)
			goto ex;
	}
	rt = 0;
ex:
	coeffmat_close(cm);
	return rt;
}

int data_coeff_matrix_csf_count(const data_id fid, size_t *n)
{
	struct coeffmat_handle cm;
	if (coeffmat_open(fid, &cm) < 0)
		return -1;

	int rt = -1;
	const htri_t has_csf =
		H5Lexists(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_CSF,
			H5P_DEFAULT);
	if (has_csf < 0)
		goto ex;
	if (!has_csf) {
		*n = 0;
		rt = 0;
		goto ex;
	}

	const hid_t cg = H5Gopen(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_CSF,
		H5P_DEFAULT);
	if (cg == H5I_INVALID_HID)
		goto ex;
	uint32_t ncomp = 0;
	if (read_u32_attr(cg, DATA_STPREP_COEFFMAT_CSF_NCOMP, &ncomp) < 0) {
		H5Gclose(cg);
		goto ex;
	}
	H5Gclose(cg);

	/*
	 * An explicit `csf/` subgroup with n_components == 0 is
	 * non-physical: the dispatcher cannot tell it apart from a
	 * single-block file with no superposition.  Reject early so
	 * the operator rebuilds the pak instead of running on a
	 * silently empty state.
	 */
	if (ncomp == 0) {
		fprintf(stderr,
			"simul.h5: /state_prep/coeff_matrix/csf present "
			"with n_components=0; remove the csf subgroup "
			"or list at least one component\n");
		rt = -EINVAL;
		goto ex;
	}

	*n = ncomp;
	rt = 0;
ex:
	coeffmat_close(cm);
	return rt;
}

int data_coeff_matrix_csf_read(const data_id fid, const size_t k,
	double *coefficient, double *C_alpha, double *C_beta)
{
	int rt = -1;

	uint32_t nqb, n_sites, n_alpha, n_beta;
	int closed_shell, tapered;
	if (data_coeff_matrix_getnums(fid, &nqb, &n_sites, &n_alpha, &n_beta,
		    &closed_shell, &tapered) < 0)
		return -1;

	if (closed_shell && C_beta != NULL)
		return -1;
	if (!closed_shell && C_beta == NULL)
		return -1;

	struct coeffmat_handle cm;
	if (coeffmat_open(fid, &cm) < 0)
		return -1;

	const hid_t cg = H5Gopen(cm.coeffmat_grp_id, DATA_STPREP_COEFFMAT_CSF,
		H5P_DEFAULT);
	if (cg == H5I_INVALID_HID)
		goto ex_csf;

	char kname[32];
	snprintf(kname, sizeof kname, "%zu", k);
	const hid_t cg_k = H5Gopen(cg, kname, H5P_DEFAULT);
	if (cg_k == H5I_INVALID_HID)
		goto ex_kgrp;

	if (read_double_attr(cg_k, DATA_STPREP_COEFFMAT_CSF_CF, coefficient)
		< 0)
		goto ex_read;
	if (read_C_dset(cg_k, DATA_STPREP_COEFFMAT_CA, n_sites, n_alpha,
		    C_alpha) < 0)
		goto ex_read;
	if (!closed_shell) {
		if (read_C_dset(cg_k, DATA_STPREP_COEFFMAT_CB, n_sites,
			    n_beta, C_beta) < 0)
			goto ex_read;
	}
	rt = 0;
ex_read:
	H5Gclose(cg_k);
ex_kgrp:
	H5Gclose(cg);
ex_csf:
	coeffmat_close(cm);
	return rt;
}
