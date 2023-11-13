#include "hdf5.h"

#include "sdat.h"
#include "log/log.h"

sdat_pauli_hamil *
sdat_pauli_hamil_create() {
    sdat_pauli_hamil *ph = malloc(sizeof(sdat_pauli_hamil));
    if (ph == NULL) {
        return NULL;
    }
    ph->coeffs = NULL;
    ph->paulis = NULL;

    return ph;
}

void
sdat_pauli_hamil_free(sdat_pauli_hamil *ph) {
    free(ph->coeffs);
    free(ph->paulis);
}


sdat_result
sdat_pauli_hamil_read(hid_t obj_id, sdat_pauli_hamil *ph) {
    sdat_result res;

    hid_t dset_coeffs_id = H5Dopen2(obj_id,
                                    SDAT_PAULI_HAMIL_COEFFS,
                                    H5P_DEFAULT);
    if (dset_coeffs_id == H5I_INVALID_HID) {
        return SDAT_ERR;
    }

    hid_t dspace_coeffs_id = H5Dget_space(dset_coeffs_id);
    hsize_t dspace_coeffs_dims[1];
    H5Sget_simple_extent_dims(dspace_coeffs_id, dspace_coeffs_dims, NULL);
    size_t num_sum_terms = dspace_coeffs_dims[0];

    double *coeffs = malloc(sizeof(double) * num_sum_terms);
    if (coeffs == NULL) {
        res = SDAT_ERR;
    }
    if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, coeffs) < 0) {
        free(coeffs);
        res = SDAT_ERR;
    }
    H5Sclose(dspace_coeffs_id);
    H5Dclose(dset_coeffs_id);

    hid_t dset_paulis_id = H5Dopen2(obj_id, SDAT_PAULI_HAMIL_PAULIS,
                                    H5P_DEFAULT);
    if (dset_paulis_id == H5I_INVALID_HID) {
        free(coeffs);
        return SDAT_ERR;
    }
    hid_t dspace_paulis_id = H5Dget_space(dset_paulis_id);
    hsize_t dspace_paulis_dims[2];
    H5Sget_simple_extent_dims(dspace_paulis_id, dspace_paulis_dims, NULL);
    if (dspace_paulis_dims[0] != num_sum_terms) {
        log_error("Number of coeffs and pauli_terms don't match");
        free(coeffs);
        H5Sclose(dspace_paulis_id);
        H5Dclose(dspace_paulis_id);
        return SDAT_ERR;
    }
    size_t num_qubits = dspace_paulis_dims[1];
    unsigned char *paulis = malloc(sizeof(unsigned char *) *
                                   num_sum_terms * num_qubits);
    H5Dread(dset_paulis_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
            H5P_DEFAULT,
            paulis);
    H5Sclose(dspace_paulis_id);
    H5Dclose(dset_paulis_id);

    ph->num_qubits = num_qubits;
    ph->num_sum_terms = num_sum_terms;
    ph->coeffs = coeffs;
    ph->paulis = paulis;

    return SDAT_OK;
}


sdat_result
sdat_time_series_read_times(hid_t grp_id,
                            sdat_time_series *ts) {

    hid_t dset_times_id = H5Dopen2(grp_id, SDAT_TIME_SERIES_TIMES,
                                   H5P_DEFAULT);
    if (dset_times_id == H5I_INVALID_HID) {
        log_error("Cannot open data set " SDAT_TIME_SERIES_TIMES);
        return SDAT_ERR;
    }
    hid_t dspace_times_id = H5Dget_space(dset_times_id);
    hsize_t dspace_times_dims[1];
    H5Sget_simple_extent_dims(dspace_times_id, dspace_times_dims, NULL);
    size_t num_steps = dspace_times_dims[0];

    double *times = malloc(sizeof(double) * num_steps);
    if (times == NULL) {
        log_fatal("Cannot allocate memory for times data set");
        H5Sclose(dspace_times_id);
        H5Dclose(dset_times_id);
        return SDAT_ERR;
    }
    H5Dread(dset_times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            times);
    H5Sclose(dspace_times_id);
    H5Dclose(dset_times_id);

    ts->num_steps = num_steps;
    ts->times = times;

    return SDAT_OK;
}


