#include "hdf5.h"

#include "data.h"

/* Group, dataset names */
#define DATA2_STATE_PREP "state_prep"
#define DATA2_STATE_PREP_MULTIDET "multidet"
#define DATA2_STATE_PREP_MULTIDET_COEFFS "coeffs"
#define DATA2_STATE_PREP_MULTIDET_DETS "dets"

/* Open, close data file */
data_id
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
		return DATA_INVALID_FID;
	}

	return file_id;
}

void
data2_close(const data_id fid)
{
	H5Fclose(fid);
}

/* --- State prep --- */
static int
state_prep_open(data_id fid, hid_t *grpid)
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
multidet_open(data_id fid, struct multidet_handle *md)
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
data2_multidet_getnums(data_id fid, size_t *num_qubits, size_t *num_dets)
{
	struct multidet_handle md;

	if (multidet_open(fid, &md) < 0)
		return -1;

	const hid_t dset_dets_id = H5Dopen2(
		md.multidet_grpid, DATA2_STATE_PREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_dets_id == H5I_INVALID_HID) {
		goto error;
	}
	const hid_t dspace_dets_id = H5Dget_space(dset_dets_id);
	hsize_t	    dspace_dets_dims[2];
	H5Sget_simple_extent_dims(dspace_dets_id, dspace_dets_dims, NULL);

	*num_dets   = dspace_dets_dims[0];
	*num_qubits = dspace_dets_dims[1];

	H5Sclose(dspace_dets_id);
	H5Dclose(dset_dets_id);
	multidet_close(md);
	return 0;

error:
	multidet_close(md);
	return -1;
}

int
multidet_read_data(
	data_id fid, _Complex double *coeffs_buf, unsigned char *dets_buf)
{
	int rc;

	struct multidet_handle md;
	if (multidet_open(fid, &md) < 0) {
		rc = -1;
		goto multidet_open_fail;
	}

	const hid_t dset_coeffs_id = H5Dopen2(md.multidet_grpid,
		DATA2_STATE_PREP_MULTIDET_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID) {
		rc = -1;
		goto coeffs_open_fail;
	}
	/* _Complex double has the same representation as double[2] */
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs_buf) < 0) {
		rc = -1;
		goto coeffs_read_fail;
	}

	const hid_t dset_dets_id = H5Dopen2(
		md.multidet_grpid, DATA2_STATE_PREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_dets_id == H5I_INVALID_HID) {
		rc = -1;
		goto dets_open_fail;
	}
	if (H5Dread(dset_dets_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, dets_buf) < 0) {
		rc = -1;
		goto dets_read_fail;
	}

	rc = 0;
dets_read_fail:
	H5Dclose(dset_dets_id);
dets_open_fail:
coeffs_read_fail:
	H5Dclose(dset_coeffs_id);
coeffs_open_fail:
	multidet_close(md);
multidet_open_fail:

	return rc;
}

int
data2_multidet_foreach(
	data_id fid, int (*op)(_Complex double, size_t, void *), void *op_data)
{
	int    rc = 0;
	size_t num_qubits, num_dets;

	if (data2_multidet_getnums(fid, &num_qubits, &num_dets) < 0)
		return -1;

	_Complex double *coeffs_buf = malloc(sizeof *coeffs_buf * num_dets);
	unsigned char	*dets_buf =
		malloc(sizeof *dets_buf * num_dets * num_qubits);
	if (!(coeffs_buf && dets_buf))
		goto error;

	/* read the content of the data file */
	if (multidet_read_data(fid, coeffs_buf, dets_buf) < 0)
		goto error;

	for (size_t i = 0; i < num_dets; i++) {
		size_t idx = 0;
		for (size_t j = 0; j < num_qubits; j++) {
			idx += dets_buf[i * num_qubits + j] << j;
		}
		rc = op(coeffs_buf[i], idx, op_data);
		if (rc != 0)
			goto exit;
	}

	goto exit;
error:
	rc = -1;
exit:
	free(coeffs_buf);
	free(dets_buf);
	return rc;
}

/* ---------------------------------------------------------------------------
 * This API is deprecated.
 */

data_id
data_file_open(const char *filename)
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
		return DATA_INVALID_FID;
	}

	return file_id;
}

void
data_file_close(const data_id fid)
{
	H5Fclose(fid);
}

void
data_state_prep_multidet_init(struct data_state_prep_multidet *dat)
{
	dat->num_qubits = 0;
	dat->num_terms	= 0;
	dat->coeffs	= NULL;
	dat->dets	= NULL;
}

void
data_state_prep_multidet_destroy(struct data_state_prep_multidet *dat)
{
	if (dat->dets) {
		free(dat->dets);
		dat->dets = NULL;
	}
	if (dat->coeffs) {
		free(dat->coeffs);
		dat->coeffs = NULL;
	}
	dat->num_terms	= 0;
	dat->num_qubits = 0;
}

int
data_state_prep_multidet_parse(
	struct data_state_prep_multidet *dat, const data_id obj_id)
{
	int res = 0;

	const hid_t dset_coeffs_id =
		H5Dopen2(obj_id, DATA_STATE_PREP_MULTIDET_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID) {
		res = -1;
		goto dset_coeffs_fail;
	}
	const hid_t dspace_coeffs_id = H5Dget_space(dset_coeffs_id);
	hsize_t	    dspace_coeffs_dims[2];
	H5Sget_simple_extent_dims(dspace_coeffs_id, dspace_coeffs_dims, NULL);
	const size_t	 num_terms = dspace_coeffs_dims[0];
	_Complex double *coeffs	   = malloc(sizeof(*coeffs) * num_terms);
	if (!coeffs) {
		res = -1;
		goto coeffs_alloc_fail;
	}
	/* _Complex double has the same representation as double[2] */
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0) {
		free(coeffs);
		coeffs = NULL;
		res    = -1;
	}

	const hid_t dset_dets_id =
		H5Dopen2(obj_id, DATA_STATE_PREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_dets_id == H5I_INVALID_HID) {
		res = -1;
		goto dset_dets_fail;
	}
	const hid_t dspace_dets_id = H5Dget_space(dset_dets_id);
	hsize_t	    dspace_dets_dims[2];
	H5Sget_simple_extent_dims(dspace_dets_id, dspace_dets_dims, NULL);
	/* It must be that:
	 * dspace_dets_dims[0] = dspace_coeffs_dims[0] = num_terms */
	const size_t   num_qubits = dspace_dets_dims[1];
	unsigned char *dets = malloc(sizeof(*dets) * num_terms * num_qubits);
	if (!dets) {
		res = -1;
		goto dets_alloc_fail;
	}
	if (H5Dread(dset_dets_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, dets) < 0) {
		free(dets);
		dets = NULL;
		res  = -1;
	}

	dat->num_qubits = num_qubits;
	dat->num_terms	= num_terms;

	dat->dets = dets;
dets_alloc_fail:
	H5Sclose(dspace_dets_id);
	H5Dclose(dset_dets_id);
dset_dets_fail:
	dat->coeffs = coeffs;
coeffs_alloc_fail:
	H5Sclose(dspace_coeffs_id);
	H5Dclose(dset_coeffs_id);
dset_coeffs_fail:
	return res;
}

void
data_state_prep_init(struct data_state_prep *dat)
{
	(void)dat;
}

void
data_state_prep_destroy(struct data_state_prep *dat)
{
	data_state_prep_multidet_destroy(&dat->multidet);
}

int
data_state_prep_parse(struct data_state_prep *dat, const data_id obj_id)
{
	int res;

	const hid_t multidet_id =
		H5Gopen2(obj_id, DATA_STATE_PREP_MULTIDET, H5P_DEFAULT);
	if (multidet_id == H5I_INVALID_HID) {
		res = -1;
		goto multidet_fail;
	}
	struct data_state_prep_multidet multidet;
	data_state_prep_multidet_init(&multidet);
	res = data_state_prep_multidet_parse(&multidet, multidet_id);
	if (res != 0) {
		data_state_prep_multidet_destroy(&multidet);
	}

	dat->multidet = multidet;
	H5Gclose(multidet_id);
multidet_fail:
	return res;
}

void
data_pauli_hamil_init(struct data_pauli_hamil *dat)
{
	dat->num_qubits = 0;
	dat->num_terms	= 0;
	dat->coeffs	= NULL;
	dat->paulis	= NULL;
}

void
data_pauli_hamil_destroy(struct data_pauli_hamil *dat)
{
	if (dat->paulis) {
		free(dat->paulis);
		dat->paulis = NULL;
	}
	if (dat->coeffs) {
		free(dat->coeffs);
		dat->coeffs = NULL;
	}
	dat->num_terms	= 0;
	dat->num_qubits = 0;
}

int
data_pauli_hamil_parse(struct data_pauli_hamil *dat, const data_id obj_id)
{
	int res = 0;

	const hid_t dset_coeffs_id =
		H5Dopen2(obj_id, DATA_PAULI_HAMIL_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID) {
		res = -1;
		goto coeffs_fail;
	}
	const hid_t dspace_coeffs_id = H5Dget_space(dset_coeffs_id);
	hsize_t	    dspace_coeffs_dims[1];
	H5Sget_simple_extent_dims(dspace_coeffs_id, dspace_coeffs_dims, NULL);
	const size_t num_terms = dspace_coeffs_dims[0];
	double	    *coeffs    = malloc(sizeof(double) * num_terms);
	if (!coeffs) {
		res = -1;
		goto coeffs_alloc_fail;
	}
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0) {
		free(coeffs);
		coeffs = NULL;
		res    = -1;
	}

	const hid_t dset_paulis_id =
		H5Dopen2(obj_id, DATA_PAULI_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_paulis_id == H5I_INVALID_HID) {
		res = -1;
		goto paulis_fail;
	}
	const hid_t dspace_paulis_id = H5Dget_space(dset_paulis_id);
	hsize_t	    dspace_paulis_dims[2];
	H5Sget_simple_extent_dims(dspace_paulis_id, dspace_paulis_dims, NULL);

	const size_t   num_qubits = dspace_paulis_dims[1];
	unsigned char *paulis =
		malloc(sizeof(unsigned char *) * num_terms * num_qubits);
	if (!paulis) {
		res = -1;
		goto paulis_alloc_fail;
	}
	H5Dread(dset_paulis_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		paulis);

	const hid_t attr_norm_id =
		H5Aopen(obj_id, DATA_PAULI_HAMIL_NORM, H5P_DEFAULT);
	if (attr_norm_id == H5I_INVALID_HID) {
		res = -1;
		goto attr_norm_fail;
	}
	double norm;
	H5Aread(attr_norm_id, H5T_NATIVE_DOUBLE, &norm);

	dat->num_qubits = num_qubits;
	dat->num_terms	= num_terms;

	dat->norm = norm;
	H5Aclose(attr_norm_id);
attr_norm_fail:
	dat->paulis = paulis;
paulis_alloc_fail:
	H5Sclose(dspace_paulis_id);
	H5Dclose(dset_paulis_id);
paulis_fail:
	dat->coeffs = coeffs;
coeffs_alloc_fail:
	H5Sclose(dspace_coeffs_id);
	H5Dclose(dset_coeffs_id);
coeffs_fail:
	return res;
}

void
data_time_series_init(struct data_time_series *dat)
{
	dat->num_steps = 0;
	dat->times     = NULL;
	dat->values    = NULL;
}

void
data_time_series_destroy(struct data_time_series *dat)
{
	if (dat->values) {
		free(dat->values);
		dat->values = NULL;
	}
	if (dat->times) {
		free(dat->times);
		dat->times = NULL;
	}
	dat->num_steps = 0;
}

int
data_time_series_parse(struct data_time_series *dat, const data_id obj_id)
{
	int res = 0;

	const hid_t dset_times_id =
		H5Dopen2(obj_id, DATA_TIME_SERIES_TIMES, H5P_DEFAULT);
	if (dset_times_id == H5I_INVALID_HID) {
		res = -1;
		goto times_fail;
	}
	const hid_t dspace_times_id = H5Dget_space(dset_times_id);
	hsize_t	    dspace_times_dims[1];
	H5Sget_simple_extent_dims(dspace_times_id, dspace_times_dims, NULL);
	const size_t num_steps = dspace_times_dims[0];
	double	    *times     = malloc(sizeof(double) * num_steps);
	if (times == NULL) {
		res = -1;
		goto times_alloc_fail;
	}
	H5Dread(dset_times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		times);

	const hid_t dset_values_id =
		H5Dopen2(obj_id, DATA_TIME_SERIES_VALUES, H5P_DEFAULT);
	if (dset_values_id == H5I_INVALID_HID) {
		res = -1;
		goto values_fail;
	}
	const hid_t dspace_values_id = H5Dget_space(dset_values_id);
	hsize_t	    dspace_values_dims[2];
	H5Sget_simple_extent_dims(dspace_values_id, dspace_values_dims, NULL);
	_Complex double *values = malloc(sizeof(*values) * num_steps);
	if (values == NULL) {
		res = -1;
		goto values_alloc_fail;
	}
	/* _Complex double has the same representation as double[2] */
	H5Dread(dset_values_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, values);

	dat->num_steps = num_steps;
	dat->values    = values;
values_alloc_fail:
	H5Sclose(dspace_values_id);
	H5Dclose(dset_values_id);
values_fail:
	dat->times = times;
times_alloc_fail:
	H5Sclose(dspace_times_id);
	H5Dclose(dset_times_id);
times_fail:
	return res;
}

int
data_time_series_write(const hid_t fid, const struct data_time_series *dat)
{
	int res = 0;

	const hid_t grp_id = H5Gopen2(fid, DATA_TIME_SERIES, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID) {
		res = -1;
		goto grp_fail;
	}
	const hid_t dset_times_id =
		H5Dopen2(grp_id, DATA_TIME_SERIES_TIMES, H5P_DEFAULT);
	if (dset_times_id == H5I_INVALID_HID) {
		res = -1;
		goto dset_times_fail;
	}
	if (H5Dwrite(dset_times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, dat->times) < 0) {
		res = -1;
	}
	H5Dclose(dset_times_id);

	const hid_t dset_values_id =
		H5Dopen2(grp_id, DATA_TIME_SERIES_VALUES, H5P_DEFAULT);
	if (dset_values_id == H5I_INVALID_HID) {
		res = -1;
		goto dset_values_fail;
	}
	/* _Complex double has the same representation of double[2] */
	if (H5Dwrite(dset_values_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, dat->values) < 0) {
		res = -1;
	}
	H5Dclose(dset_values_id);

dset_values_fail:
dset_times_fail:
	H5Gclose(grp_id);
grp_fail:
	return res;
}

void
data_init(struct data *dat)
{
	(void)dat;
}

void
data_destroy(struct data *dat)
{
	data_state_prep_destroy(&dat->state_prep);
	data_pauli_hamil_destroy(&dat->pauli_hamil);
	data_time_series_destroy(&dat->time_series);
}

int
data_parse(struct data *dat, const data_id fid)
{
	int res;

	const hid_t state_prep_id = H5Gopen2(fid, DATA_STATE_PREP, H5P_DEFAULT);
	if (state_prep_id == H5I_INVALID_HID) {
		res = -1;
		goto state_prep_fail;
	}
	struct data_state_prep state_prep;
	data_state_prep_init(&state_prep);
	res = data_state_prep_parse(&state_prep, state_prep_id);
	if (res != 0) {
		data_state_prep_destroy(&state_prep);
	}

	const hid_t pauli_hamil_id =
		H5Gopen2(fid, DATA_PAULI_HAMIL, H5P_DEFAULT);
	if (pauli_hamil_id == H5I_INVALID_HID) {
		res = -1;
		goto pauli_hamil_fail;
	}
	struct data_pauli_hamil pauli_hamil;
	data_pauli_hamil_init(&pauli_hamil);
	res = data_pauli_hamil_parse(&pauli_hamil, pauli_hamil_id);
	if (res != 0) {
		data_pauli_hamil_destroy(&pauli_hamil);
	}

	const hid_t time_series_id =
		H5Gopen2(fid, DATA_TIME_SERIES, H5P_DEFAULT);
	if (time_series_id == H5I_INVALID_HID) {
		res = -1;
		goto time_series_fail;
	}
	struct data_time_series time_series;
	data_time_series_init(&time_series);
	res = data_time_series_parse(&time_series, time_series_id);
	if (res != 0) {
		data_time_series_destroy(&time_series);
	}

	dat->time_series = time_series;
	H5Gclose(time_series_id);
time_series_fail:
	dat->pauli_hamil = pauli_hamil;
	H5Gclose(pauli_hamil_id);
pauli_hamil_fail:
	dat->state_prep = state_prep;
	H5Gclose(state_prep_id);
state_prep_fail:
	return res;
}

/* --------------------------------------------------------------------------
 */
