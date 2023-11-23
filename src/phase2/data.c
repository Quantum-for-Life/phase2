#include "hdf5.h"

#include "data.h"
#include "log.h"

void data_pauli_hamil_init(struct data_pauli_hamil* dat)
{
        dat->num_qubits = 0;
        dat->num_sum_terms = 0;
        dat->coeffs = NULL;
        dat->paulis = NULL;
}

void data_pauli_hamil_destroy(struct data_pauli_hamil* dat)
{
        free(dat->coeffs);
        dat->coeffs = NULL;
        free(dat->paulis);
        dat->paulis = NULL;
}

int data_pauli_hamil_read(struct data_pauli_hamil* dat, const hid_t obj_id)
{
        int res = DATA_OK;

        const hid_t grp_id = H5Gopen2(obj_id, DATA_PAULI_HAMIL, H5P_DEFAULT);
        if (grp_id == H5I_INVALID_HID)
        {
                res = DATA_ERR;
                goto grp_fail;
        }
        const hid_t dset_coeffs_id = H5Dopen2(grp_id,
                                              DATA_PAULI_HAMIL_COEFFS,
                                              H5P_DEFAULT);
        if (dset_coeffs_id == H5I_INVALID_HID)
        {
                res = DATA_ERR;
                goto dset_coeffs_fail;
        }
        const hid_t dspace_coeffs_id = H5Dget_space(dset_coeffs_id);
        hsize_t dspace_coeffs_dims[1];
        H5Sget_simple_extent_dims(dspace_coeffs_id, dspace_coeffs_dims, NULL);
        dat->num_sum_terms = dspace_coeffs_dims[0];
        double* coeffs = malloc(sizeof(double) * dat->num_sum_terms);
        if (coeffs == NULL)
        {
                res = DATA_ERR;
                goto coeffs_fail;
        }
        dat->coeffs = coeffs;
        if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE,
                    H5S_ALL, H5S_ALL, H5P_DEFAULT, coeffs) < 0)
        {
                free(coeffs);
                res = DATA_ERR;
        }

        const hid_t dset_paulis_id = H5Dopen2(grp_id,
                                              DATA_PAULI_HAMIL_PAULIS,
                                              H5P_DEFAULT);
        if (dset_paulis_id == H5I_INVALID_HID)
        {
                res = DATA_ERR;
                goto dset_paulis_fail;
        }
        const hid_t dspace_paulis_id = H5Dget_space(dset_paulis_id);
        hsize_t dspace_paulis_dims[2];
        H5Sget_simple_extent_dims(dspace_paulis_id, dspace_paulis_dims, NULL);
        if (dspace_paulis_dims[0] != dat->num_sum_terms)
        {
                res = DATA_ERR;
                goto dim_mismatch;
        }
        dat->num_qubits = dspace_paulis_dims[1];
        unsigned char* paulis = malloc(sizeof(unsigned char*) *
                dat->num_sum_terms * dat->num_qubits);
        if (paulis == NULL)
        {
                res = DATA_ERR;
                goto paulis_fail;
        }
        H5Dread(dset_paulis_id,
                H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
                H5P_DEFAULT,
                paulis);
        dat->paulis = paulis;

paulis_fail:
dim_mismatch:
        H5Sclose(dspace_paulis_id);
        H5Dclose(dset_paulis_id);
dset_paulis_fail:
coeffs_fail:
        H5Sclose(dspace_coeffs_id);
        H5Dclose(dset_coeffs_id);
dset_coeffs_fail:
        H5Gclose(grp_id);
grp_fail:

        return res;
}


void data_time_series_init(struct data_time_series* dat)
{
        dat->num_steps = 0;
        dat->times = NULL;
        dat->values = NULL;
}

void data_time_series_destroy(struct data_time_series* dat)
{
        free(dat->times);
        dat->times = NULL;
        free(dat->values);
        dat->values = NULL;
}

int data_time_series_read(struct data_time_series* dat, const hid_t obj_id)
{
        int res = DATA_OK;

        const hid_t grp_id = H5Gopen2(obj_id, DATA_TIME_SERIES, H5P_DEFAULT);
        if (grp_id == H5I_INVALID_HID)
        {
                res = DATA_ERR;
                goto grp_fail;
        }
        const hid_t dset_times_id = H5Dopen2(grp_id, DATA_TIME_SERIES_TIMES,
                                             H5P_DEFAULT);
        if (dset_times_id == H5I_INVALID_HID)
        {
                res = DATA_ERR;
                goto dset_times_fail;
        }
        const hid_t dspace_times_id = H5Dget_space(dset_times_id);
        hsize_t dspace_times_dims[1];
        H5Sget_simple_extent_dims(dspace_times_id, dspace_times_dims, NULL);
        dat->num_steps = dspace_times_dims[0];
        double* times = malloc(sizeof(double) * dat->num_steps);
        if (times == NULL)
        {
                res = DATA_ERR;
                goto times_fail;
        }
        H5Dread(dset_times_id,
                H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                times);
        dat->times = times;

        const hid_t dset_values_id = H5Dopen2(grp_id,
                                              DATA_TIME_SERIES_VALUES,
                                              H5P_DEFAULT);
        if (dset_values_id == H5I_INVALID_HID)
        {
                res = DATA_ERR;
                goto dset_values_fail;
        }
        const hid_t dspace_values_id = H5Dget_space(dset_values_id);
        hsize_t dspace_values_dims[2];
        H5Sget_simple_extent_dims(dspace_values_id, dspace_values_dims, NULL);
        if (dspace_values_dims[0] != dat->num_steps ||
                dspace_values_dims[1] != 2)
        {
                res = DATA_ERR;
                goto values_dims_mismatch;
        }
        double* values = malloc(sizeof(double) * dat->num_steps * 2);
        if (values == NULL)
        {
                res = DATA_ERR;
                goto values_fail;
        }
        H5Dread(dset_values_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                H5P_DEFAULT,
                values);
        dat->values = values;

values_fail:
values_dims_mismatch:
        H5Sclose(dspace_values_id);
        H5Dclose(dset_values_id);
dset_values_fail:
times_fail:
        H5Sclose(dspace_times_id);
        H5Dclose(dset_times_id);
dset_times_fail:
        H5Gclose(grp_id);
grp_fail:

        return res;
}

int data_time_series_write(const struct data_time_series dat,
                           const hid_t obj_id)
{
        int res = DATA_OK;

        const hid_t grp_id = H5Gopen2(obj_id, DATA_TIME_SERIES, H5P_DEFAULT);
        if (grp_id == H5I_INVALID_HID)
        {
                res = DATA_ERR;
                goto grp_fail;
        }
        const hid_t dset_times_id = H5Dopen2(grp_id,
                                             DATA_TIME_SERIES_TIMES,
                                             H5P_DEFAULT);
        if (dset_times_id == H5I_INVALID_HID)
        {
                res = DATA_ERR;
                goto dset_times_fail;
        }
        if (H5Dwrite(dset_times_id,
                     H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, dat.times) < 0)
        {
                res = DATA_ERR;
        }
        H5Dclose(dset_times_id);

        const hid_t dset_values_id = H5Dopen2(grp_id,
                                              DATA_TIME_SERIES_VALUES,
                                              H5P_DEFAULT);
        if (dset_values_id == H5I_INVALID_HID)
        {
                res = DATA_ERR;
                goto dset_values_fail;
        }
        if (H5Dwrite(dset_values_id,
                     H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, dat.values) < 0)
        {
                res = DATA_ERR;
        }
        H5Dclose(dset_values_id);

dset_values_fail:
dset_times_fail:
        H5Gclose(grp_id);
grp_fail:

        return res;
}
