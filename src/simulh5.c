#include "hdf5.h"

#include "simulh5.h"
#include "log/log.h"


typedef struct {
    size_t num_qubits;
    size_t num_sum_terms;
    double *coeffs;
    unsigned char *paulis;
} simulh5_grp_pauli_hamil;

typedef struct {
    size_t num_steps;
    double *times;
    double *values_real;
    double *values_imag;
} simulh5_grp_time_series;


simulh5_result
simulh5_grp_pauli_hamil_read(hid_t grp_id, simulh5_grp_pauli_hamil *ph) {

    log_debug("Read group " SIMULH5_GRP_PAULI_HAMIL);
    hid_t dset_coeffs_id = H5Dopen2(grp_id,
                                    SIMULH5_GRP_PAULI_HAMIL_COEFFS,
                                    H5P_DEFAULT);
    if (dset_coeffs_id == H5I_INVALID_HID) {
        return SIMULH5_ERR;
    }

    hid_t dspace_coeffs_id = H5Dget_space(dset_coeffs_id);
    hsize_t dspace_coeffs_dims[1];
    H5Sget_simple_extent_dims(dspace_coeffs_id, dspace_coeffs_dims, NULL);
    size_t num_sum_terms = dspace_coeffs_dims[0];
    log_debug("Number of coefficients: %zu", num_sum_terms);

    double *coeffs = malloc(sizeof(double) * num_sum_terms);
    if (coeffs == NULL) {
        return SIMULH5_ERR;
    }
    if (H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, coeffs) < 0) {
        free(coeffs);
        return SIMULH5_ERR;
    }
    H5Sclose(dspace_coeffs_id);
    H5Dclose(dset_coeffs_id);

    hid_t dset_paulis_id = H5Dopen2(grp_id, SIMULH5_GRP_PAULI_HAMIL_PAULIS,
                                    H5P_DEFAULT);
    if (dset_paulis_id == H5I_INVALID_HID) {
        free(coeffs);
        return SIMULH5_ERR;
    }
    hid_t dspace_paulis_id = H5Dget_space(dset_paulis_id);
    hsize_t dspace_paulis_dims[2];
    H5Sget_simple_extent_dims(dspace_paulis_id, dspace_paulis_dims, NULL);
    if (dspace_paulis_dims[0] != num_sum_terms) {
        log_error("Number of coeffs and pauli_terms don't match");
        free(coeffs);
        H5Sclose(dspace_paulis_id);
        H5Dclose(dspace_paulis_id);
        return SIMULH5_ERR;
    }
    size_t num_qubits = dspace_paulis_dims[1];
    log_debug("Number of qubits %zu", num_qubits);
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

    return SIMULH5_OK;
}


simulh5_result
simulh5_grp_time_series_read_times(hid_t grp_id,
                                   simulh5_grp_time_series *ts) {

    hid_t dset_times_id = H5Dopen2(grp_id, SIMULH5_GRP_TIME_SERIES_TIMES,
                                   H5P_DEFAULT);
    if (dset_times_id == H5I_INVALID_HID) {
        log_error("Cannot open data set " SIMULH5_GRP_TIME_SERIES_TIMES);
        return SIMULH5_ERR;
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
        return SIMULH5_ERR;
    }
    H5Dread(dset_times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            times);
    H5Sclose(dspace_times_id);
    H5Dclose(dset_times_id);

    ts->num_steps = num_steps;
    ts->times = times;

    return SIMULH5_OK;
}


simulh5 *simulh5_create() {
    simulh5 *sh = malloc(sizeof(simulh5));
    if (sh == NULL) {
        return NULL;
    }

    sh->time_series.values_real = NULL;
    sh->time_series.values_imag = NULL;

    return sh;
}

simulh5_result
simulh5_read(simulh5 *sh, hid_t obj_id) {
    simulh5_result res;

    hid_t pauli_hamil_id = H5Gopen2(obj_id, SIMULH5_GRP_PAULI_HAMIL,
                                    H5P_DEFAULT);
    simulh5_grp_pauli_hamil pauli_hamil;
    res = simulh5_grp_pauli_hamil_read(pauli_hamil_id, &pauli_hamil);
    H5Gclose(pauli_hamil_id);
    if (res != SIMULH5_OK) {
        simulh5_free(sh);
        return res;
    }

    hid_t time_series_id = H5Gopen2(obj_id, SIMULH5_GRP_TIME_SERIES,
                                    H5P_DEFAULT);
    simulh5_grp_time_series time_series;
    res = simulh5_grp_time_series_read_times(time_series_id, &time_series);
    H5Gclose(time_series_id);
    if (res != SIMULH5_OK) {
        simulh5_free(sh);
        return res;
    }

    sh->pauli_hamil.num_qubits = pauli_hamil.num_qubits;
    sh->pauli_hamil.num_sum_terms = pauli_hamil.num_sum_terms;
    sh->pauli_hamil.coeffs = pauli_hamil.coeffs;
    sh->pauli_hamil.paulis = pauli_hamil.paulis;

    sh->time_series.num_steps = time_series.num_steps;
    sh->time_series.times = time_series.times;
    sh->time_series.values_real = NULL;
    sh->time_series.values_imag = NULL;

    return SIMULH5_OK;
}


void simulh5_free(simulh5 *sh) {
    free(sh->pauli_hamil.coeffs);
    free(sh->pauli_hamil.paulis);
    free(sh->time_series.times);
    free(sh->time_series.values_real);
    free(sh->time_series.values_imag);
    free(sh);
}
