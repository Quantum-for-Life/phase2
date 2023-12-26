#include "hdf5.h"

#include "data.h"

/* Group, dataset names */
#define DATA2_STATE_PREP "state_prep"
#define DATA2_STATE_PREP_MULTIDET "multidet"
#define DATA2_STATE_PREP_MULTIDET_COEFFS "coeffs"
#define DATA2_STATE_PREP_MULTIDET_DETS "dets"

#define DATA2_PAULI_HAMIL "pauli_hamil"
#define DATA2_PAULI_HAMIL_COEFFS "coeffs"
#define DATA2_PAULI_HAMIL_PAULIS "paulis"
#define DATA2_PAULI_HAMIL_NORM "normalization"

#define DATA2_TIME_SERIES "time_series"
#define DATA2_TIME_SERIES_TIMES "times"
#define DATA2_TIME_SERIES_VALUES "values"

/* Open, close data file */
data2_id
data2_open(const char *filename)
{
	hid_t file_id, access_plist;
#ifdef DISTRIBUTED
	// init MPI environment
	int initialized;
	MPI_Initialized(&initialized);
	if (!initialized) {
		MPI_Init(NULL, NULL);
	}
	access_plist = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(access_plist, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
	access_plist = H5P_DEFAULT;
#endif
	file_id = H5Fopen(filename, H5F_ACC_RDWR, access_plist);
	if (file_id == H5I_INVALID_HID) {
		return DATA2_INVALID_FID;
	}

	return file_id;
}

void
data2_close(const data2_id fid)
{
	H5Fclose(fid);
}

/* --- State prep --- */
static int
state_prep_open(data2_id fid, hid_t *grpid)
{
	hid_t hid = H5Gopen2(fid, DATA2_STATE_PREP, H5P_DEFAULT);
	if (hid == H5I_INVALID_HID)
		return -1;
	*grpid = hid;

	return 0;
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
multidet_open(data2_id fid, struct multidet_handle *md)
{
	hid_t sp_id, md_id;

	if (state_prep_open(fid, &sp_id) < 0)
		return -1;
	md_id = H5Gopen2(sp_id, DATA2_STATE_PREP_MULTIDET, H5P_DEFAULT);
	if (md_id == H5I_INVALID_HID) {
		state_prep_close(sp_id);
		return -1;
	}
	md->state_prep_grpid = sp_id;
	md->multidet_grpid   = md_id;

	return 0;
}

static void
multidet_close(struct multidet_handle md)
{
	H5Gclose(md.multidet_grpid);
	state_prep_close(md.state_prep_grpid);
}

int
data2_multidet_getnums(data2_id fid, size_t *num_qubits, size_t *num_dets)
{
	struct multidet_handle md;

	if (multidet_open(fid, &md) < 0)
		goto err_open;

	const hid_t dset_id = H5Dopen2(
		md.multidet_grpid, DATA2_STATE_PREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID)
		goto err_md_open;

	const hid_t dsp_id = H5Dget_space(dset_id);
	hsize_t	    dsp_dims[2];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 2)
		goto err_dims;

	*num_dets   = dsp_dims[0];
	*num_qubits = dsp_dims[1];

	H5Sclose(dsp_id);
	H5Dclose(dset_id);
	multidet_close(md);
	return 0;

err_dims:
	H5Sclose(dsp_id);
	H5Dclose(dset_id);
err_md_open:
	multidet_close(md);
err_open:
	return -1;
}

static int
multidet_read_data(data2_id fid, _Complex double *coeffs, unsigned char *dets)
{
	struct multidet_handle md;
	if (multidet_open(fid, &md) < 0)
		goto err_multidet_open;

	const hid_t dset_coeffs_id = H5Dopen2(md.multidet_grpid,
		DATA2_STATE_PREP_MULTIDET_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID)
		goto err_coeffs_open;
	/* _Complex double has the same representation as double[2] */
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0)
		goto err_coeffs_read;

	const hid_t dset_dets_id = H5Dopen2(
		md.multidet_grpid, DATA2_STATE_PREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_dets_id == H5I_INVALID_HID)
		goto err_dets_open;
	if (H5Dread(dset_dets_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, dets) < 0)
		goto err_dets_read;

	H5Dclose(dset_dets_id);
	H5Dclose(dset_coeffs_id);
	multidet_close(md);
	return 0;

err_dets_read:
	H5Dclose(dset_dets_id);
err_dets_open:
err_coeffs_read:
	H5Dclose(dset_coeffs_id);
err_coeffs_open:
	multidet_close(md);
err_multidet_open:
	return -1;
}

int
data2_multidet_foreach(
	data2_id fid, int (*op)(_Complex double, size_t, void *), void *op_data)
{
	int    rc = 0;
	size_t num_qubits, num_dets;

	if (data2_multidet_getnums(fid, &num_qubits, &num_dets) < 0)
		goto err_getnums;

	_Complex double *coeffs_buf = malloc(sizeof *coeffs_buf * num_dets);
	if (!coeffs_buf)
		goto err_alloc_coeffs;
	unsigned char *dets_buf =
		malloc(sizeof *dets_buf * num_dets * num_qubits);
	if (!dets_buf)
		goto err_alloc_dets;

	/* Read the content of the data file */
	if (multidet_read_data(fid, coeffs_buf, dets_buf) < 0)
		goto err_data_read;

	for (size_t i = 0; i < num_dets; i++) {
		size_t idx = 0;
		for (size_t j = 0; j < num_qubits; j++) {
			idx += dets_buf[i * num_qubits + j] << j;
		}
		rc = op(coeffs_buf[i], idx, op_data);
		/* This isn't an error, but rather the user telling us to
		   short-circuit the iteration. */
		if (rc != 0)
			break;
	}

	free(dets_buf);
	free(coeffs_buf);
	return rc;

err_data_read:
	free(dets_buf);
err_alloc_dets:
	free(coeffs_buf);
err_alloc_coeffs:
err_getnums:
	return -1;
}

static int
hamil_open(data2_id fid, hid_t *grpid)
{
	const hid_t hamil_id = H5Gopen2(fid, DATA2_PAULI_HAMIL, H5P_DEFAULT);
	if (hamil_id == H5I_INVALID_HID)
		return -1;
	*grpid = hamil_id;

	return 0;
}

static void
hamil_close(hid_t grpid)
{
	H5Gclose(grpid);
}

int
data2_hamil_getnums(data2_id fid, size_t *num_qubits, size_t *num_terms)
{
	hid_t grpid;
	if (hamil_open(fid, &grpid) < 0)
		goto err_open;

	const hid_t dset_id =
		H5Dopen2(grpid, DATA2_PAULI_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID)
		goto err_read;

	const hid_t dsp_id = H5Dget_space(dset_id);
	hsize_t	    dsp_dims[2];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 2)
		goto err_dims;

	*num_terms  = dsp_dims[0];
	*num_qubits = dsp_dims[1];

	H5Sclose(dsp_id);
	H5Dclose(dset_id);
	hamil_close(grpid);
	return 0;

err_dims:
	H5Sclose(dsp_id);
	H5Dclose(dset_id);
err_read:
	hamil_close(grpid);
err_open:
	return -1;
}

int
data2_hamil_getnorm(data2_id fid, double *norm)
{
	hid_t grpid;
	if (hamil_open(fid, &grpid) < 0)
		goto err_open;

	const hid_t attr_norm_id =
		H5Aopen(grpid, DATA2_PAULI_HAMIL_NORM, H5P_DEFAULT);
	if (attr_norm_id == H5I_INVALID_HID)
		goto err_attr_open;
	double n;
	if (H5Aread(attr_norm_id, H5T_NATIVE_DOUBLE, &n) < 0)
		goto err_attr_read;
	*norm = n;

	H5Aclose(attr_norm_id);
	hamil_close(grpid);
	return 0;

err_attr_read:
	H5Aclose(attr_norm_id);
err_attr_open:
	hamil_close(grpid);
err_open:
	return -1;
}

static int
hamil_read_data(data2_id fid, double *coeffs, unsigned char *paulis)
{
	hid_t grpid;
	if (hamil_open(fid, &grpid) < 0)
		goto err_hamil_open;

	const hid_t dset_coeffs_id =
		H5Dopen2(grpid, DATA2_PAULI_HAMIL_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID)
		goto err_coeffs_open;
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0)
		goto err_coeffs_read;

	const hid_t dset_paulis_id =
		H5Dopen2(grpid, DATA2_PAULI_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_paulis_id == H5I_INVALID_HID)
		goto err_paulis_open;
	if (H5Dread(dset_paulis_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, paulis) < 0)
		goto err_pauli_read;

	H5Dclose(dset_paulis_id);
	H5Dclose(dset_coeffs_id);
	hamil_close(grpid);
	return 0;

err_pauli_read:
	H5Dclose(dset_paulis_id);
err_paulis_open:
err_coeffs_read:
	H5Dclose(dset_coeffs_id);
err_coeffs_open:
	hamil_close(grpid);
err_hamil_open:
	return -1;
}

int
data2_hamil_foreach(
	data2_id fid, int (*op)(double, unsigned char *, void *), void *op_data)
{
	int	       rc = 0;
	unsigned char *paulis, *paustr;
	double	      *coeffs;
	size_t	       num_qubits, num_terms;

	if (data2_hamil_getnums(fid, &num_qubits, &num_terms) < 0)
		goto err_getnums;
	coeffs = malloc(sizeof *coeffs * num_terms);
	if (!coeffs)
		goto err_coeffs_alloc;
	paulis = malloc(sizeof *paulis * num_qubits * num_terms);
	if (!paulis)
		goto err_paulis_alloc;
	if (hamil_read_data(fid, coeffs, paulis) < 0)
		goto err_hamil_read;
	paustr = malloc(sizeof *paustr * num_qubits);
	if (!paustr)
		goto err_paustr_alloc;

	for (size_t i = 0; i < num_terms; i++) {
		for (size_t j = 0; j < num_qubits; j++) {
			paustr[j] = paulis[i * num_qubits + j];
		}
		rc = op(coeffs[i], paustr, op_data);
		if (rc != 0)
			break;
	}

	free(paustr);
	free(paulis);
	free(coeffs);
	return rc;

err_paustr_alloc:
err_hamil_read:
	free(paulis);
err_paulis_alloc:
	free(coeffs);
err_coeffs_alloc:
err_getnums:
	return -1;
}

static int
times_open(data2_id fid, hid_t *grpid)
{
	const hid_t id = H5Gopen2(fid, DATA2_TIME_SERIES, H5P_DEFAULT);
	if (id == H5I_INVALID_HID)
		return -1;

	*grpid = id;
	return 0;
}

static void
times_close(hid_t grpid)
{
	H5Gclose(grpid);
}

int
data2_times_getnums(data2_id fid, size_t *num_steps)
{
	hid_t grpid;
	if (times_open(fid, &grpid) < 0)
		goto err_open;

	const hid_t dset_id =
		H5Dopen2(grpid, DATA2_TIME_SERIES_TIMES, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID)
		goto err_dset_open;

	const hid_t dsp_id = H5Dget_space(dset_id);
	hsize_t	    dsp_dims[1];
	if (H5Sget_simple_extent_dims(dsp_id, dsp_dims, NULL) != 1)
		goto err_dims;

	*num_steps = dsp_dims[0];

	H5Sclose(dsp_id);
	H5Dclose(dset_id);
	times_close(grpid);
	return 0;

err_dims:
	H5Sclose(dsp_id);
	H5Dclose(dset_id);
err_dset_open:
	times_close(grpid);
err_open:
	return -1;
}

static int
times_read_data(data2_id fid, double *times, _Complex double *values)
{
	hid_t grpid;
	if (times_open(fid, &grpid) < 0)
		goto err_open;

	const hid_t times_id =
		H5Dopen2(grpid, DATA2_TIME_SERIES_TIMES, H5P_DEFAULT);
	if (times_id == H5I_INVALID_HID)
		goto err_times;
	if (H5Dread(times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    times) < 0)
		goto err_times_read;

	const hid_t values_id =
		H5Dopen2(grpid, DATA2_TIME_SERIES_VALUES, H5P_DEFAULT);
	if (values_id == H5I_INVALID_HID)
		goto err_values;
	/* _Complex double has the same representation of double[2] */
	if (H5Dread(values_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    values) < 0)
		goto err_values_read;

	H5Dclose(values_id);
	H5Dclose(times_id);
	times_close(grpid);
	return 0;

err_values_read:
	H5Dclose(values_id);
err_values:
err_times_read:
	H5Dclose(times_id);
err_times:
	times_close(grpid);
err_open:
	return -1;
}

static int
times_write_data(data2_id fid, double *times, _Complex double *values)
{
	hid_t grpid;
	if (times_open(fid, &grpid) < 0)
		goto err_open;

	const hid_t times_id =
		H5Dopen2(grpid, DATA2_TIME_SERIES_TIMES, H5P_DEFAULT);
	if (times_id == H5I_INVALID_HID)
		goto err_times;
	if (H5Dwrite(times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    times) < 0)
		goto err_times_write;

	const hid_t values_id =
		H5Dopen2(grpid, DATA2_TIME_SERIES_VALUES, H5P_DEFAULT);
	if (values_id == H5I_INVALID_HID)
		goto err_values;
	/* _Complex double has the same representation of double[2] */
	if (H5Dwrite(values_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, values) < 0)
		goto err_values_write;

	H5Dclose(values_id);
	H5Dclose(times_id);
	times_close(grpid);
	return 0;

err_values_write:
	H5Dclose(values_id);
err_values:
err_times_write:
	H5Dclose(times_id);
err_times:
	times_close(grpid);
err_open:
	return -1;
}

int
data2_times_foreach(
	data2_id fid, int (*op)(double, _Complex double, void *), void *op_data)
{
	int rc = 0;

	double		*times;
	_Complex double *values;
	size_t		 num_steps;

	if (data2_times_getnums(fid, &num_steps) < 0)
		goto err_getnums;
	times = malloc(sizeof *times * num_steps);
	if (!times)
		goto err_times_alloc;
	values = malloc(sizeof *values * num_steps);
	if (!values)
		goto err_values_alloc;
	if (times_read_data(fid, times, values) < 0)
		goto err_times_read;

	for (size_t i = 0; i < num_steps; i++) {
		rc = op(times[i], values[i], op_data);
		if (rc != 0)
			break;
	}

	free(values);
	free(times);
	return rc;

err_times_read:
	free(values);
err_values_alloc:
	free(times);
err_times_alloc:
err_getnums:
	return -1;
}

int
data2_times_update(data2_id fid, int (*op)(double *, _Complex double *, void *),
	void *op_data)
{
	int rc = 0;

	double		*times;
	_Complex double *values;
	size_t		 num_steps;

	if (data2_times_getnums(fid, &num_steps) < 0)
		goto err_getnums;
	times = malloc(sizeof *times * num_steps);
	if (!times)
		goto err_times_alloc;
	values = malloc(sizeof *values * num_steps);
	if (!values)
		goto err_values_alloc;
	if (times_read_data(fid, times, values) < 0)
		goto err_times_read;
	for (size_t i = 0; i < num_steps; i++) {
		rc = op(&times[i], &values[i], op_data);
		if (rc != 0)
			break;
	}
	if (times_write_data(fid, times, values) < 0)
		goto err_times_write;

	free(values);
	free(times);
	return rc;

err_times_write:
err_times_read:
	free(values);
err_values_alloc:
	free(times);
err_times_alloc:
err_getnums:
	return -1;
}
