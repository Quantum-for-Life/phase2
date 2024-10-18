#include "c23_compat.h"
#include <stddef.h>
#include <stdint.h>

#include "hdf5.h"

#include "phase2/data.h"

/* Group, dataset names */
#define DATA_STPREP "state_prep"
#define DATA_STPREP_MULTIDET "multidet"
#define DATA_STPREP_MULTIDET_COEFFS "coeffs"
#define DATA_STPREP_MULTIDET_DETS "dets"

#define DATA_HAMIL "pauli_hamil"
#define DATA_HAMIL_COEFFS "coeffs"
#define DATA_HAMIL_PAULIS "paulis"
#define DATA_HAMIL_NORM "normalization"

#define DATA_CIRCTROTT "circ_trott"
#define DATA_CIRCTROTT_TIMEFACTOR "time_factor"
#define DATA_CIRCTROTT_VALUES "values"

#define DATA_CIRCQDRIFT "circ_qdrift"
#define DATA_CIRCQDRIFT_STEPSIZE "step_size"
#define DATA_CIRCQDRIFT_NUMSAMPLES "num_samples"
#define DATA_CIRCQDRIFT_DEPTH "depth"
#define DATA_CIRCQDRIFT_VALUES "values"

/* Open, close data file */
data_id data_open(const char *filename)
{
	// init MPI environment
	int initialized;
	MPI_Initialized(&initialized);
	if (!initialized)
		MPI_Init(NULL, NULL);
	const hid_t access_plist = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(access_plist, MPI_COMM_WORLD, MPI_INFO_NULL);

	const hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, access_plist);
	if (file_id == H5I_INVALID_HID)
		return DATA_INVALID_FID;

	return file_id;
}

void data_close(const data_id fid)
{
	H5Fclose(fid);
}

static int data_group_open(const hid_t fid, hid_t *grpid, const char *name)
{
	const hid_t hid = H5Gopen2(fid, name, H5P_DEFAULT);
	if (hid == H5I_INVALID_HID)
		return -1;
	*grpid = hid;

	return 0;
}

static void data_group_close(const hid_t grpid)
{
	H5Gclose(grpid);
}

static int data_read_attr(hid_t fid, const char *name, hid_t type_id, void *buf)
{
	int rt = 0;

	hid_t attr_id = H5Aopen(fid, name, H5P_DEFAULT);
	if (attr_id == H5I_INVALID_HID)
		return -1;
	if (H5Aread(attr_id, type_id, buf) < 0)
		rt = -1;
	H5Aclose(attr_id);

	return rt;
}

struct multidet_handle {
	hid_t state_prep_grpid;
	hid_t multidet_grpid;
};

static int multidet_open(const data_id fid, struct multidet_handle *md)
{
	hid_t sp_id, md_id;

	if (data_group_open(fid, &sp_id, DATA_STPREP) < 0)
		return -1;
	if (data_group_open(sp_id, &md_id, DATA_STPREP_MULTIDET) < 0) {
		data_group_close(sp_id);
		return -1;
	}

	md->state_prep_grpid = sp_id;
	md->multidet_grpid = md_id;

	return 0;
}

static void multidet_close(struct multidet_handle md)
{
	data_group_close(md.multidet_grpid);
	data_group_close(md.state_prep_grpid);
}

int data_multidet_getnums(data_id fid, uint32_t *num_qubits, size_t *num_dets)
{
	int rt = -1;

	struct multidet_handle md;
	if (multidet_open(fid, &md) < 0)
		return -1;

	const hid_t dset_id = H5Dopen2(
		md.multidet_grpid, DATA_STPREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID)
		goto exit_md_open;

	const hid_t dsp_id = H5Dget_space(dset_id);
	hsize_t dsp_dims[2];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 2)
		goto exit_dims;

	*num_dets = dsp_dims[0];
	*num_qubits = (uint32_t)dsp_dims[1];
	rt = 0;

exit_dims:
	H5Sclose(dsp_id);
	H5Dclose(dset_id);
exit_md_open:
	multidet_close(md);
	return rt;
}

static int multidet_read_data(data_id fid, double *coeffs, unsigned char *dets)
{
	int rt = -1;

	struct multidet_handle md;
	if (multidet_open(fid, &md) < 0)
		goto exit_multidet_open;

	const hid_t dset_coeffs_id = H5Dopen2(
		md.multidet_grpid, DATA_STPREP_MULTIDET_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID)
		goto exit_coeffs_open;
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0)
		goto exit_coeffs_read;

	const hid_t dset_dets_id = H5Dopen2(
		md.multidet_grpid, DATA_STPREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_dets_id == H5I_INVALID_HID)
		goto exit_dets_open;
	if (H5Dread(dset_dets_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, dets) < 0)
		goto exit_dets_read;

	rt = 0;

exit_dets_read:
	H5Dclose(dset_dets_id);
exit_dets_open:
exit_coeffs_read:
	H5Dclose(dset_coeffs_id);
exit_coeffs_open:
	multidet_close(md);
exit_multidet_open:
	return rt;
}

int data_multidet_foreach(data_id fid,
	int (*op)(double coeff[2], uint64_t idx, void *), void *op_data)
{
	int rt = -1, rc = 0;
	uint32_t num_qubits;	
	size_t num_dets;

	if (data_multidet_getnums(fid, &num_qubits, &num_dets) < 0)
		goto exit_getnums;

	double *coeffs_buf = malloc(sizeof *coeffs_buf * 2 * num_dets);
	if (coeffs_buf == NULL)
		goto exit_alloc_coeffs;
	unsigned char *dets_buf =
		malloc(sizeof *dets_buf * num_dets * num_qubits);
	if (dets_buf == NULL)
		goto exit_alloc_dets;

	/* Read the content of the data file */
	if (multidet_read_data(fid, coeffs_buf, dets_buf) < 0)
		goto exit_data_read;

	for (size_t i = 0; i < num_dets; i++) {
		uint64_t it_idx = 0;
		for (size_t j = 0; j < num_qubits; j++) {
			it_idx += dets_buf[i * num_qubits + j] << j;
		}
		double it_coeff[2] = {

			coeffs_buf[2 * i], coeffs_buf[2 * i + 1]
		};
		rc = op(it_coeff, it_idx, op_data);
		/* This isn't an error, but rather the user telling us to
		   short-circuit the iteration. */
		if (rc != 0)
			break;
	}
	rt = rc;

exit_data_read:
	free(dets_buf);
exit_alloc_dets:
	free(coeffs_buf);
exit_alloc_coeffs:
exit_getnums:
	return rt;
}

int data_hamil_getnums(data_id fid, uint32_t *num_qubits, size_t *num_terms)
{
	int rt = -1;

	hid_t grpid;
	if (data_group_open(fid, &grpid, DATA_HAMIL) < 0)
		goto exit_open;

	const hid_t dset_id = H5Dopen2(grpid, DATA_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID)
		goto exit_read;

	const hid_t dsp_id = H5Dget_space(dset_id);
	hsize_t dsp_dims[2];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 2)
		goto exit_dims;

	*num_terms = dsp_dims[0];
	*num_qubits = (uint32_t)dsp_dims[1];
	rt = 0;

exit_dims:
	H5Sclose(dsp_id);
	H5Dclose(dset_id);
exit_read:
	data_group_close(grpid);
exit_open:
	return rt;
}

int data_hamil_getnorm(data_id fid, double *norm)
{
	int rt = -1;

	hid_t grpid;
	if (data_group_open(fid, &grpid, DATA_HAMIL) < 0)
		goto exit_open;

	const hid_t attr_norm_id = H5Aopen(grpid, DATA_HAMIL_NORM, H5P_DEFAULT);
	if (attr_norm_id == H5I_INVALID_HID)
		goto exit_attr_open;

	double n;
	if (H5Aread(attr_norm_id, H5T_NATIVE_DOUBLE, &n) < 0)
		goto exit_attr_read;

	*norm = n;
	rt = 0;

exit_attr_read:
	H5Aclose(attr_norm_id);
exit_attr_open:
	data_group_close(grpid);
exit_open:
	return rt;
}

static int hamil_read_data(data_id fid, double *coeffs, unsigned char *paulis)
{
	int rt = -1;

	hid_t grpid;
	if (data_group_open(fid, &grpid, DATA_HAMIL) < 0)
		return -1;

	const hid_t dset_coeffs_id =
		H5Dopen2(grpid, DATA_HAMIL_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID)
		goto exit_coeffs_open;

	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0)
		goto exit_coeffs_read;

	const hid_t dset_paulis_id =
		H5Dopen2(grpid, DATA_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_paulis_id == H5I_INVALID_HID)
		goto exit_paulis_open;

	if (H5Dread(dset_paulis_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, paulis) < 0)
		goto exit_pauli_read;

	rt = 0;

exit_pauli_read:
	H5Dclose(dset_paulis_id);
exit_paulis_open:
exit_coeffs_read:
	H5Dclose(dset_coeffs_id);
exit_coeffs_open:
	data_group_close(grpid);

	return rt;
}

int data_hamil_foreach(const data_id fid,
	int (*op)(double, unsigned char *, void *), void *op_data)
{
	int rt = -1, rc = 0;
	uint32_t num_qubits;
	size_t num_terms;

	if (data_hamil_getnums(fid, &num_qubits, &num_terms) < 0)
		return -1;

	double *coeffs = malloc(sizeof *coeffs * num_terms);
	if (coeffs == NULL)
		goto exit_coeffs_alloc;

	unsigned char *paulis = malloc(sizeof *paulis * num_qubits * num_terms);
	if (paulis == NULL)
		goto exit_paulis_alloc;
	if (hamil_read_data(fid, coeffs, paulis) < 0)
		goto exit_hamil_read;
	unsigned char *paustr = malloc(sizeof *paustr * num_qubits);
	if (paustr == NULL)
		goto exit_paustr_alloc;

	for (size_t i = 0; i < num_terms; i++) {
		for (size_t j = 0; j < num_qubits; j++)
			paustr[j] = paulis[i * num_qubits + j];
		rc = op(coeffs[i], paustr, op_data);
		if (rc != 0)
			break;
	}
	free(paustr);
	rt = rc;

exit_paustr_alloc:
exit_hamil_read:
	free(paulis);
exit_paulis_alloc:
	free(coeffs);
exit_coeffs_alloc:
	return rt;
}

int data_circ_trott_getttrs(data_id fid, double *delta)
{
	int rt = 0;

	hid_t grpid;
	if (data_group_open(fid, &grpid, DATA_CIRCTROTT) < 0)
		return -1;

	rt += data_read_attr(
		grpid, DATA_CIRCTROTT_TIMEFACTOR, H5T_NATIVE_DOUBLE, delta);

	data_group_close(grpid);

	return rt;
}

static int data_circ_write_values(const char *grp_name, const char *dset_name,
	const data_id fid, double *values[2], const size_t num_values)
{
	int rt = -1;

	hid_t grpid, dspace, dset;
	double *val_cont = malloc(sizeof(double) * 2 * num_values);
	if (val_cont == NULL)
		return -1;
	for (size_t i = 0; i < num_values; i++) {
		val_cont[2 * i] = values[0][i];
		val_cont[2 * i + 1] = values[1][i];
	}

	if (data_group_open(fid, &grpid, grp_name) < 0)
		goto exit_open;
	if ((dspace = H5Screate_simple(
		     2, (hsize_t[]){ num_values, 2 }, NULL)) == H5I_INVALID_HID)
		goto exit_fspace;

	if ((dset = H5Dcreate2(grpid, dset_name, H5T_IEEE_F64LE, dspace,
		     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) == H5I_INVALID_HID)
		goto exit_dset;

	if (H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, dspace, H5P_DEFAULT,
		    val_cont) < 0)
		goto exit_dset_write;

	rt = 0;

exit_dset_write:
	H5Dclose(dset);
exit_dset:
	H5Sclose(dspace);
exit_fspace:
	data_group_close(grpid);
exit_open:
	free(val_cont);

	return rt;
}

int data_circ_trott_write_values(
	const data_id fid, double *values[2], const size_t num_values)
{
	return data_circ_write_values(
		DATA_CIRCTROTT, DATA_CIRCTROTT_VALUES, fid, values, num_values);
}

int data_circ_trott_read_values_test(
	const data_id fid, double *values[2], const size_t num_values)
{
	int rt = -1;

	hid_t grpid;
	double *val_cont = malloc(sizeof(double) * 2 * num_values);
	if (val_cont == NULL)
		return -1;

	if (data_group_open(fid, &grpid, DATA_CIRCTROTT) < 0)
		goto exit_open;
	const hid_t dset = H5Dopen2(grpid, DATA_CIRCTROTT_VALUES, H5P_DEFAULT);
	if (dset == H5I_INVALID_HID)
		goto exit_dset;
	if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    val_cont) < 0)
		goto exit_read;

	for (size_t i = 0; i < num_values; i++) {
		values[0][i] = val_cont[2 * i];
		values[1][i] = val_cont[2 * i + 1];
	}

	rt = 0;

exit_read:
	H5Dclose(dset);
exit_dset:
	data_group_close(grpid);
exit_open:
	free(val_cont);

	return rt;
}

int data_circ_qdrift_getattrs(const data_id fid, size_t *num_samples,
	double *step_size, size_t *depth)
{
	int rt = 0;

	hid_t grpid;
	if (data_group_open(fid, &grpid, DATA_CIRCQDRIFT) < 0)
		return -1;

	rt += data_read_attr(grpid, DATA_CIRCQDRIFT_NUMSAMPLES,
		H5T_NATIVE_UINT64, num_samples);
	rt += data_read_attr(
		grpid, DATA_CIRCQDRIFT_STEPSIZE, H5T_NATIVE_DOUBLE, step_size);
	rt += data_read_attr(
		grpid, DATA_CIRCQDRIFT_DEPTH, H5T_NATIVE_UINT64, depth);

	data_group_close(grpid);

	return rt;
}

int data_circ_qdrift_write_values(
	const data_id fid, double *values[2], const size_t num_values)
{
	return data_circ_write_values(DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_VALUES,
		fid, values, num_values);
}
