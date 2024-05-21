#include <stdint.h>

#include "hdf5.h"

#include "data.h"

enum {
	DATA_OK = 0,
	DATA_EDIM,
	DATA_EFILE,
	DATA_ENOMEM,
	DATA_ESTPREP,
	DATA_EMULDET,
	DATA_EHAMIL,
	DATA_ETROTT,
};

/* Group, dataset names */
#define DATA_STATE_PREP "state_prep"
#define DATA_STATE_PREP_MULTIDET "multidet"
#define DATA_STATE_PREP_MULTIDET_COEFFS "coeffs"
#define DATA_STATE_PREP_MULTIDET_DETS "dets"

#define DATA_PAULI_HAMIL "pauli_hamil"
#define DATA_PAULI_HAMIL_COEFFS "coeffs"
#define DATA_PAULI_HAMIL_PAULIS "paulis"
#define DATA_PAULI_HAMIL_NORM "normalization"

#define DATA_CIRC_QDRIFT "circ_qdrift"
#define DATA_CIRC_QDRIFT_STEP_SIZE "step_size"
#define DATA_CIRC_QDRIFT_NUM_SAMPLES "num_samples"
#define DATA_CIRC_QDRIFT_DEPTH "depth"
#define DATA_CIRC_QDRIFT_VALUES "values"

/* Open, close data file */
data_id
data_open(const char *filename)
{
	// init MPI environment
	int initialized;
	MPI_Initialized(&initialized);
	if (!initialized) {
		MPI_Init(NULL, NULL);
	}
	const hid_t access_plist = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(access_plist, MPI_COMM_WORLD, MPI_INFO_NULL);

	const hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, access_plist);
	if (file_id == H5I_INVALID_HID) {
		return DATA_INVALID_FID;
	}

	return file_id;
}

void
data_close(const data_id fid)
{
	H5Fclose(fid);
}

/* --- State prep --- */
static int
state_prep_open(const data_id fid, hid_t *grpid)
{
	const hid_t hid = H5Gopen2(fid, DATA_STATE_PREP, H5P_DEFAULT);
	if (hid == H5I_INVALID_HID)
		return -DATA_ESTPREP;
	*grpid = hid;

	return DATA_OK;
}

static void
state_prep_close(hid_t grpid)
{
	H5Gclose(grpid);
}

/* --- Multidet --- */

struct multidet_handle {
	hid_t state_prep_grpid;
	hid_t multidet_grpid;
};

static int
multidet_open(const data_id fid, struct multidet_handle *md)
{
	hid_t sp_id;

	if (state_prep_open(fid, &sp_id) < 0)
		return -DATA_ESTPREP;
	const hid_t md_id =
		H5Gopen2(sp_id, DATA_STATE_PREP_MULTIDET, H5P_DEFAULT);
	if (md_id == H5I_INVALID_HID) {
		state_prep_close(sp_id);
		return -DATA_EMULDET;
	}
	md->state_prep_grpid = sp_id;
	md->multidet_grpid   = md_id;

	return DATA_OK;
}

static void
multidet_close(struct multidet_handle md)
{
	H5Gclose(md.multidet_grpid);
	state_prep_close(md.state_prep_grpid);
}

int
data_multidet_getnums(data_id fid, size_t *num_qubits, size_t *num_dets)
{
	int ret = DATA_OK;

	struct multidet_handle md;

	if (multidet_open(fid, &md) < 0) {
		ret = -DATA_EFILE;
		goto err_open;
	}

	const hid_t dset_id = H5Dopen2(
		md.multidet_grpid, DATA_STATE_PREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID) {
		ret = -DATA_EMULDET;
		goto err_md_open;
	}

	const hid_t dsp_id = H5Dget_space(dset_id);
	hsize_t	    dsp_dims[2];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 2) {
		ret = -DATA_EDIM;
		goto err_dims;
	}

	*num_dets   = dsp_dims[0];
	*num_qubits = dsp_dims[1];

err_dims:
	H5Sclose(dsp_id);
	H5Dclose(dset_id);
err_md_open:
	multidet_close(md);
err_open:
	return ret;
}

static int
multidet_read_data(data_id fid, double *coeffs, unsigned char *dets)
{
	int ret = DATA_OK;

	struct multidet_handle md;
	if (multidet_open(fid, &md) < 0) {
		ret = -DATA_EMULDET;
		goto err_multidet_open;
	}

	const hid_t dset_coeffs_id = H5Dopen2(md.multidet_grpid,
		DATA_STATE_PREP_MULTIDET_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID) {
		ret = DATA_EMULDET;
		goto err_coeffs_open;
	}
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0) {
		ret = -DATA_EMULDET;
		goto err_coeffs_read;
	}

	const hid_t dset_dets_id = H5Dopen2(
		md.multidet_grpid, DATA_STATE_PREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_dets_id == H5I_INVALID_HID) {
		ret = -DATA_EMULDET;
		goto err_dets_open;
	}
	if (H5Dread(dset_dets_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, dets) < 0) {
		ret = -DATA_EMULDET;
		goto err_dets_read;
	}

err_dets_read:
	H5Dclose(dset_dets_id);
err_dets_open:
err_coeffs_read:
	H5Dclose(dset_coeffs_id);
err_coeffs_open:
	multidet_close(md);
err_multidet_open:

	return ret;
}

int
data_multidet_foreach(data_id fid,
	int (*op)(double coeff[2], uint64_t idx, void *), void *op_data)
{
	int    ret, rc = 0;
	size_t num_qubits, num_dets;

	if (data_multidet_getnums(fid, &num_qubits, &num_dets) < 0) {
		ret = -DATA_EMULDET;
		goto err_getnums;
	}

	double *coeffs_buf = malloc(sizeof *coeffs_buf * 2 * num_dets);
	if (coeffs_buf == NULL) {
		ret = -DATA_ENOMEM;
		goto err_alloc_coeffs;
	}
	unsigned char *dets_buf =
		malloc(sizeof *dets_buf * num_dets * num_qubits);
	if (dets_buf == NULL) {
		ret = -DATA_ENOMEM;
		goto err_alloc_dets;
	}

	/* Read the content of the data file */
	if (multidet_read_data(fid, coeffs_buf, dets_buf) < 0) {
		ret = -DATA_EMULDET;
		goto err_data_read;
	}

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
	ret = rc;

err_data_read:
	free(dets_buf);
err_alloc_dets:
	free(coeffs_buf);
err_alloc_coeffs:
err_getnums:

	return ret;
}

static int
hamil_open(data_id fid, hid_t *grpid)
{
	const hid_t hamil_id = H5Gopen2(fid, DATA_PAULI_HAMIL, H5P_DEFAULT);
	if (hamil_id == H5I_INVALID_HID)
		return -DATA_EHAMIL;

	*grpid = hamil_id;

	return DATA_OK;
}

static void
hamil_close(hid_t grpid)
{
	H5Gclose(grpid);
}

int
data_hamil_getnums(data_id fid, size_t *num_qubits, size_t *num_terms)
{
	int ret = DATA_OK;

	hid_t grpid;
	if (hamil_open(fid, &grpid) < 0) {
		ret = -DATA_EHAMIL;
		goto err_open;
	}

	const hid_t dset_id =
		H5Dopen2(grpid, DATA_PAULI_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID) {
		ret = -DATA_EHAMIL;
		goto err_read;
	}

	const hid_t dsp_id = H5Dget_space(dset_id);
	hsize_t	    dsp_dims[2];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 2) {
		ret = -DATA_EDIM;
		goto err_dims;
	}

	*num_terms  = dsp_dims[0];
	*num_qubits = dsp_dims[1];

err_dims:
	H5Sclose(dsp_id);
	H5Dclose(dset_id);
err_read:
	hamil_close(grpid);
err_open:

	return ret;
}

int
data_hamil_getnorm(data_id fid, double *norm)
{
	int ret = DATA_OK;

	hid_t grpid;
	if (hamil_open(fid, &grpid) < 0) {
		ret = -DATA_EHAMIL;
		goto err_open;
	}

	const hid_t attr_norm_id =
		H5Aopen(grpid, DATA_PAULI_HAMIL_NORM, H5P_DEFAULT);
	if (attr_norm_id == H5I_INVALID_HID) {
		ret = -DATA_EHAMIL;
		goto err_attr_open;
	}
	double n;
	if (H5Aread(attr_norm_id, H5T_NATIVE_DOUBLE, &n) < 0) {
		ret = -DATA_EHAMIL;
		goto err_attr_read;
	}
	*norm = n;

err_attr_read:
	H5Aclose(attr_norm_id);
err_attr_open:
	hamil_close(grpid);
err_open:

	return ret;
}

static int
hamil_read_data(data_id fid, double *coeffs, unsigned char *paulis)
{
	int ret = DATA_OK;

	hid_t grpid;
	if (hamil_open(fid, &grpid) < 0)
		goto err_hamil_open;

	const hid_t dset_coeffs_id =
		H5Dopen2(grpid, DATA_PAULI_HAMIL_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID) {
		ret = -DATA_EHAMIL;
		goto err_coeffs_open;
	}
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0) {
		ret = -DATA_EHAMIL;
		goto err_coeffs_read;
	}

	const hid_t dset_paulis_id =
		H5Dopen2(grpid, DATA_PAULI_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_paulis_id == H5I_INVALID_HID) {
		ret = -DATA_EHAMIL;
		goto err_paulis_open;
	}
	if (H5Dread(dset_paulis_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, paulis) < 0) {
		ret = -DATA_EHAMIL;
		goto err_pauli_read;
	}

err_pauli_read:
	H5Dclose(dset_paulis_id);
err_paulis_open:
err_coeffs_read:
	H5Dclose(dset_coeffs_id);
err_coeffs_open:
	hamil_close(grpid);
err_hamil_open:

	return ret;
}

int
data_hamil_foreach(const data_id fid,
	int (*op)(double, unsigned char *, void *), void *op_data)
{
	int    ret, rc = 0;
	size_t num_qubits, num_terms;

	if (data_hamil_getnums(fid, &num_qubits, &num_terms) < 0) {
		ret = -DATA_EHAMIL;
		goto err_getnums;
	}
	double *coeffs = malloc(sizeof *coeffs * num_terms);
	if (coeffs == NULL) {
		ret = -DATA_EHAMIL;
		goto err_coeffs_alloc;
	}
	unsigned char *paulis = malloc(sizeof *paulis * num_qubits * num_terms);
	if (paulis == NULL) {
		ret = -DATA_EHAMIL;
		goto err_paulis_alloc;
	}
	if (hamil_read_data(fid, coeffs, paulis) < 0) {
		ret = -DATA_EHAMIL;
		goto err_hamil_read;
	}
	unsigned char *paustr = malloc(sizeof *paustr * num_qubits);
	if (paustr == NULL) {
		ret = -DATA_ENOMEM;
		goto err_paustr_alloc;
	}

	for (size_t i = 0; i < num_terms; i++) {
		for (size_t j = 0; j < num_qubits; j++) {
			paustr[j] = paulis[i * num_qubits + j];
		}
		rc = op(coeffs[i], paustr, op_data);
		if (rc != 0)
			break;
	}
	free(paustr);
	ret = rc;

err_paustr_alloc:
err_hamil_read:
	free(paulis);
err_paulis_alloc:
	free(coeffs);
err_coeffs_alloc:
err_getnums:
	return ret;
}

static int
trotter_open(data_id fid, hid_t *grpid)
{
	const hid_t id = H5Gopen2(fid, DATA_CIRC_QDRIFT, H5P_DEFAULT);
	if (id == H5I_INVALID_HID)
		return -DATA_ETROTT;

	*grpid = id;
	return DATA_OK;
}

static void
trotter_close(hid_t grpid)
{
	H5Gclose(grpid);
}

int
data_circ_qdrift_get_factor(data_id fid, double *step_size)
{
	int ret = DATA_OK;

	hid_t grpid;
	if (trotter_open(fid, &grpid) < 0) {
		ret = -DATA_ETROTT;
		goto err_open;
	}

	const hid_t attr_norm_id =
		H5Aopen(grpid, DATA_CIRC_QDRIFT_STEP_SIZE, H5P_DEFAULT);
	if (attr_norm_id == H5I_INVALID_HID) {
		ret = -DATA_ETROTT;
		goto err_attr_open;
	}
	double n;
	if (H5Aread(attr_norm_id, H5T_NATIVE_DOUBLE, &n) < 0) {
		ret = -DATA_ETROTT;
		goto err_attr_read;
	}
	*step_size = n;

err_attr_read:
	H5Aclose(attr_norm_id);
err_attr_open:
	trotter_close(grpid);
err_open:

	return ret;
}


int
data_circ_qdrift_get_num_samples(data_id fid, size_t *num_samples)
{
	int ret = DATA_OK;

	hid_t grpid;
	if (trotter_open(fid, &grpid) < 0) {
		ret = -DATA_ETROTT;
		goto err_open;
	}

	const hid_t attr_norm_id =
		H5Aopen(grpid, DATA_CIRC_QDRIFT_NUM_SAMPLES, H5P_DEFAULT);
	if (attr_norm_id == H5I_INVALID_HID) {
		ret = -DATA_ETROTT;
		goto err_attr_open;
	}
	uint64_t n;
	if (H5Aread(attr_norm_id, H5T_NATIVE_UINT64, &n) < 0) {
		ret = -DATA_ETROTT;
		goto err_attr_read;
	}
	*num_samples = (size_t)n;

err_attr_read:
	H5Aclose(attr_norm_id);
err_attr_open:
	trotter_close(grpid);
err_open:

	return ret;
}

int
data_circ_qdrift_get_depth(data_id fid, size_t *depth)
{
	int ret = DATA_OK;

	hid_t grpid;
	if (trotter_open(fid, &grpid) < 0) {
		ret = -DATA_ETROTT;
		goto err_open;
	}

	const hid_t attr_norm_id =
		H5Aopen(grpid, DATA_CIRC_QDRIFT_DEPTH, H5P_DEFAULT);
	if (attr_norm_id == H5I_INVALID_HID) {
		ret = -DATA_ETROTT;
		goto err_attr_open;
	}
	uint64_t n;
	if (H5Aread(attr_norm_id, H5T_NATIVE_UINT64, &n) < 0) {
		ret = -DATA_ETROTT;
		goto err_attr_read;
	}
	*depth = (size_t)n;

err_attr_read:
	H5Aclose(attr_norm_id);
err_attr_open:
	trotter_close(grpid);
err_open:

	return ret;
}

int
data_circ_qdrift_write_values(data_id fid, double *values[2], size_t num_values)
{
	int ret = DATA_OK;

	hid_t grpid, dspace, dset;

	double *val_cont = malloc(sizeof(double) * 2 * num_values);
	if (val_cont == NULL)
		return -1;
	for (size_t i = 0; i < num_values; i++) {
		val_cont[2 * i]	    = values[0][i];
		val_cont[2 * i + 1] = values[1][i];
	}

	if (trotter_open(fid, &grpid) < 0) {
		ret = -DATA_ETROTT;
		goto err_open;
	}
	if ((dspace = H5Screate_simple(2, (hsize_t[]){ num_values, 2 },
		     NULL)) == H5I_INVALID_HID) {
		goto err_fspace;
	}
	if ((dset = H5Dcreate2(grpid, DATA_CIRC_QDRIFT_VALUES, H5T_IEEE_F64LE,
		     dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) ==
		H5I_INVALID_HID) {
		ret = -DATA_ETROTT;
		goto err_dset;
	}
	if (H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, dspace, H5P_DEFAULT,
		    val_cont) < 0) {
		ret = -DATA_ETROTT;
		goto err_dset_write;
	}

err_dset_write:
	H5Dclose(dset);
err_dset:
	H5Sclose(dspace);
err_fspace:
	trotter_close(grpid);
err_open:
	free(val_cont);

	return ret;
}

int
data_trotter_read_values_test(data_id fid, double *values[2], size_t num_values)
{
	int ret = DATA_OK;

	hid_t grpid;

	double *val_cont = malloc(sizeof(double) * 2 * num_values);
	if (val_cont == NULL)
		return -DATA_ENOMEM;

	if (trotter_open(fid, &grpid) < 0) {
		ret = -DATA_ETROTT;
		goto err_open;
	}
	const hid_t dset =
		H5Dopen2(grpid, DATA_CIRC_QDRIFT_VALUES, H5P_DEFAULT);
	if (dset == H5I_INVALID_HID) {
		ret = -DATA_ETROTT;
		goto err_dset;
	}
	if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    val_cont) < 0) {
		ret = -DATA_ETROTT;
		goto err_read;
	}

	for (size_t i = 0; i < num_values; i++) {
		values[0][i] = val_cont[2 * i];
		values[1][i] = val_cont[2 * i + 1];
	}

err_read:
	H5Dclose(dset);
err_dset:
	trotter_close(grpid);
err_open:
	free(val_cont);
	return ret;
}
