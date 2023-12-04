#include "hdf5.h"

#ifdef DISTRIBUTED

#include "log.h"

#endif

#include "data.h"

data_id data_file_open(const char *filename)
{
	hid_t file_id, access_plist;
#ifdef DISTRIBUTED
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

void data_file_close(data_id fid)
{
	H5Fclose(fid);
}

void data_state_prep_multidet_init(struct data_state_prep_multidet *dat)
{
	dat->num_qubits = 0;
	dat->num_terms = 0;
	dat->coeffs = NULL;
	dat->dets = NULL;
}

void data_state_prep_multidet_destroy(struct data_state_prep_multidet *dat)
{
	if (dat->dets) {
		free(dat->dets);
		dat->dets = NULL;
	}
	if (dat->coeffs) {
		free(dat->coeffs);
		dat->coeffs = NULL;
	}
	dat->num_terms = 0;
	dat->num_qubits = 0;
}

int data_state_prep_multidet_parse(struct data_state_prep_multidet *dat,
				   data_id obj_id)
{
	int res = DATA_OK;

	const hid_t dset_coeffs_id =
		H5Dopen2(obj_id, DATA_STATE_PREP_MULTIDET_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto dset_coeffs_fail;
	}
	const hid_t dspace_coeffs_id = H5Dget_space(dset_coeffs_id);
	hsize_t dspace_coeffs_dims[2];
	H5Sget_simple_extent_dims(dspace_coeffs_id, dspace_coeffs_dims, NULL);
	const size_t num_terms = dspace_coeffs_dims[0];
	double *coeffs = malloc(sizeof(*coeffs) * num_terms * 2);
	if (!coeffs) {
		res = DATA_ERR;
		goto coeffs_alloc_fail;
	}
	if (H5Dread(dset_coeffs_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0) {
		free(coeffs);
		coeffs = NULL;
		res = DATA_ERR;
	}

	const hid_t dset_dets_id =
		H5Dopen2(obj_id, DATA_STATE_PREP_MULTIDET_DETS, H5P_DEFAULT);
	if (dset_dets_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto dset_dets_fail;
	}
	const hid_t dspace_dets_id = H5Dget_space(dset_dets_id);
	hsize_t dspace_dets_dims[2];
	H5Sget_simple_extent_dims(dspace_dets_id, dspace_dets_dims, NULL);
	/* It must be that:
         * dspace_dets_dims[0] = dspace_coeffs_dims[0] = num_terms */
	const size_t num_qubits = dspace_dets_dims[1];
	unsigned char *dets = malloc(sizeof(*dets) * num_terms * num_qubits);
	if (!dets) {
		res = DATA_ERR;
		goto dets_alloc_fail;
	}
	if (H5Dread(dset_dets_id, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    dets) < 0) {
		free(dets);
		dets = NULL;
		res = DATA_ERR;
	}

	dat->num_qubits = num_qubits;
	dat->num_terms = num_terms;

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

void data_state_prep_init(struct data_state_prep *dat)
{
	(void)dat;
}

void data_state_prep_destroy(struct data_state_prep *dat)
{
	data_state_prep_multidet_destroy(&dat->multidet);
}

int data_state_prep_parse(struct data_state_prep *dat, data_id obj_id)
{
	int res;

	const hid_t multidet_id =
		H5Gopen2(obj_id, DATA_STATE_PREP_MULTIDET, H5P_DEFAULT);
	if (multidet_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto multidet_fail;
	}
	struct data_state_prep_multidet multidet;
	data_state_prep_multidet_init(&multidet);
	res = data_state_prep_multidet_parse(&multidet, multidet_id);
	if (res != DATA_OK) {
		data_state_prep_multidet_destroy(&multidet);
	}

	dat->multidet = multidet;
	H5Gclose(multidet_id);
multidet_fail:
	return res;
}

void data_pauli_hamil_init(struct data_pauli_hamil *dat)
{
	dat->num_qubits = 0;
	dat->num_terms = 0;
	dat->coeffs = NULL;
	dat->paulis = NULL;
}

void data_pauli_hamil_destroy(struct data_pauli_hamil *dat)
{
	if (dat->paulis) {
		free(dat->paulis);
		dat->paulis = NULL;
	}
	if (dat->coeffs) {
		free(dat->coeffs);
		dat->coeffs = NULL;
	}
	dat->num_terms = 0;
	dat->num_qubits = 0;
}

int data_pauli_hamil_parse(struct data_pauli_hamil *dat, data_id obj_id)
{
	int res = DATA_OK;

	const hid_t dset_coeffs_id =
		H5Dopen2(obj_id, DATA_PAULI_HAMIL_COEFFS, H5P_DEFAULT);
	if (dset_coeffs_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto coeffs_fail;
	}
	const hid_t dspace_coeffs_id = H5Dget_space(dset_coeffs_id);
	hsize_t dspace_coeffs_dims[1];
	H5Sget_simple_extent_dims(dspace_coeffs_id, dspace_coeffs_dims, NULL);
	const size_t num_terms = dspace_coeffs_dims[0];
	double *coeffs = malloc(sizeof(double) * num_terms);
	if (!coeffs) {
		res = DATA_ERR;
		goto coeffs_alloc_fail;
	}
	if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, coeffs) < 0) {
		free(coeffs);
		coeffs = NULL;
		res = DATA_ERR;
	}

	const hid_t dset_paulis_id =
		H5Dopen2(obj_id, DATA_PAULI_HAMIL_PAULIS, H5P_DEFAULT);
	if (dset_paulis_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto paulis_fail;
	}
	const hid_t dspace_paulis_id = H5Dget_space(dset_paulis_id);
	hsize_t dspace_paulis_dims[2];
	H5Sget_simple_extent_dims(dspace_paulis_id, dspace_paulis_dims, NULL);

	const size_t num_qubits = dspace_paulis_dims[1];
	unsigned char *paulis =
		malloc(sizeof(unsigned char *) * num_terms * num_qubits);
	if (!paulis) {
		res = DATA_ERR;
		goto paulis_alloc_fail;
	}
	H5Dread(dset_paulis_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		paulis);

	const hid_t attr_norm_id =
		H5Aopen(obj_id, DATA_PAULI_HAMIL_NORM, H5P_DEFAULT);
	if (attr_norm_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto attr_norm_fail;
	}
	double norm;
	H5Aread(attr_norm_id, H5T_IEEE_F64LE, &norm);

	dat->num_qubits = num_qubits;
	dat->num_terms = num_terms;

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

void data_time_series_init(struct data_time_series *dat)
{
	dat->num_steps = 0;
	dat->times = NULL;
	dat->values = NULL;
}

void data_time_series_destroy(struct data_time_series *dat)
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

int data_time_series_parse(struct data_time_series *dat, data_id obj_id)
{
	int res = DATA_OK;

	const hid_t dset_times_id =
		H5Dopen2(obj_id, DATA_TIME_SERIES_TIMES, H5P_DEFAULT);
	if (dset_times_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto times_fail;
	}
	const hid_t dspace_times_id = H5Dget_space(dset_times_id);
	hsize_t dspace_times_dims[1];
	H5Sget_simple_extent_dims(dspace_times_id, dspace_times_dims, NULL);
	const size_t num_steps = dspace_times_dims[0];
	double *times = malloc(sizeof(double) * num_steps);
	if (times == NULL) {
		res = DATA_ERR;
		goto times_alloc_fail;
	}
	H5Dread(dset_times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		times);

	const hid_t dset_values_id =
		H5Dopen2(obj_id, DATA_TIME_SERIES_VALUES, H5P_DEFAULT);
	if (dset_values_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto values_fail;
	}
	const hid_t dspace_values_id = H5Dget_space(dset_values_id);
	hsize_t dspace_values_dims[2];
	H5Sget_simple_extent_dims(dspace_values_id, dspace_values_dims, NULL);
	double *values = malloc(sizeof(double) * num_steps * 2);
	if (values == NULL) {
		res = DATA_ERR;
		goto values_alloc_fail;
	}
	H5Dread(dset_values_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, values);

	dat->num_steps = num_steps;

	dat->values = values;
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

int data_time_series_write(const hid_t fid, const struct data_time_series *dat)
{
	int res = DATA_OK;

	const hid_t grp_id = H5Gopen2(fid, DATA_TIME_SERIES, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto grp_fail;
	}
	const hid_t dset_times_id =
		H5Dopen2(grp_id, DATA_TIME_SERIES_TIMES, H5P_DEFAULT);
	if (dset_times_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto dset_times_fail;
	}
	if (H5Dwrite(dset_times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, dat->times) < 0) {
		res = DATA_ERR;
	}
	H5Dclose(dset_times_id);

	const hid_t dset_values_id =
		H5Dopen2(grp_id, DATA_TIME_SERIES_VALUES, H5P_DEFAULT);
	if (dset_values_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto dset_values_fail;
	}
	if (H5Dwrite(dset_values_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, dat->values) < 0) {
		res = DATA_ERR;
	}
	H5Dclose(dset_values_id);

dset_values_fail:
dset_times_fail:
	H5Gclose(grp_id);
grp_fail:
	return res;
}

void data_init(struct data *dat)
{
	(void)dat;
}

void data_destroy(struct data *dat)
{
	data_state_prep_destroy(&dat->state_prep);
	data_pauli_hamil_destroy(&dat->pauli_hamil);
	data_time_series_destroy(&dat->time_series);
}

int data_parse(struct data *dat, data_id fid)
{
	int res;

	const hid_t state_prep_id = H5Gopen2(fid, DATA_STATE_PREP, H5P_DEFAULT);
	if (state_prep_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto state_prep_fail;
	}
	struct data_state_prep state_prep;
	data_state_prep_init(&state_prep);
	res = data_state_prep_parse(&state_prep, state_prep_id);
	if (res != DATA_OK) {
		data_state_prep_destroy(&state_prep);
	}

	const hid_t pauli_hamil_id =
		H5Gopen2(fid, DATA_PAULI_HAMIL, H5P_DEFAULT);
	if (pauli_hamil_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto pauli_hamil_fail;
	}
	struct data_pauli_hamil pauli_hamil;
	data_pauli_hamil_init(&pauli_hamil);
	res = data_pauli_hamil_parse(&pauli_hamil, pauli_hamil_id);
	if (res != DATA_OK) {
		data_pauli_hamil_destroy(&pauli_hamil);
	}

	const hid_t time_series_id =
		H5Gopen2(fid, DATA_TIME_SERIES, H5P_DEFAULT);
	if (time_series_id == H5I_INVALID_HID) {
		res = DATA_ERR;
		goto time_series_fail;
	}
	struct data_time_series time_series;
	data_time_series_init(&time_series);
	res = data_time_series_parse(&time_series, time_series_id);
	if (res != DATA_OK) {
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
