//
// Created by mm on 11/12/23.
//

#include "hdf5.h"

#include "h5.h"
#include "log/log.h"

int
h5_read_attr(hid_t obj_id,
             const char *attr_name,
             hid_t type_id,
             void *value) {
    __label__ fail_attr;
    int ret_val = EXIT_SUCCESS;

    hid_t attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT);
    if (attr_id == H5I_INVALID_HID) {
        ret_val = EXIT_FAILURE;
        goto fail_attr;
    }
    if (H5Aread(attr_id, type_id, value) < 0) {
        ret_val = EXIT_FAILURE;
    }

    H5Aclose(attr_id);
    fail_attr:
    return ret_val;
}


int
h5_read_dset(hid_t obj_id, const char *dset_name, hid_t type_id, void *value) {

    __label__ fail_dset;
    int ret_val = EXIT_SUCCESS;

    hid_t dset_id = H5Dopen2(obj_id, dset_name, H5P_DEFAULT);
    if (dset_id == H5I_INVALID_HID) {
        ret_val = EXIT_FAILURE;
        goto fail_dset;
    }
    if (H5Dread(dset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, value) < 0) {
        ret_val = EXIT_FAILURE;
    }
    H5Dclose(dset_id);
    fail_dset:
    return ret_val;
}

int
h5_grp_hamiltonian_read(hid_t group_id, h5_grp_hamiltonian *gh) {

    int ret_val = h5_read_attr(group_id, "num_qubits", H5T_NATIVE_INT,
                               &gh->num_qubits);
    if (ret_val != EXIT_SUCCESS) {
        return ret_val;
    }
    ret_val = h5_read_attr(group_id, "num_sum_terms", H5T_NATIVE_INT,
                           &gh->num_sum_terms);
    if (ret_val != EXIT_SUCCESS) {
        return ret_val;
    }
    double *coeffs = malloc(sizeof(double) * gh->num_qubits);
    if (coeffs == NULL) {
        return EXIT_FAILURE;
    }
    ret_val = h5_read_dset(group_id, "coeffs", H5T_NATIVE_DOUBLE, coeffs);
    if (ret_val != EXIT_SUCCESS) {
        free(coeffs);
        return ret_val;
    }
    int *paulis = malloc(sizeof(int) * gh->num_qubits *
                         gh->num_sum_terms);
    if (paulis == NULL) {
        free(coeffs);
        return EXIT_FAILURE;
    }
    ret_val = h5_read_dset(group_id, "paulis", H5T_NATIVE_INT, paulis);

    if (ret_val != EXIT_SUCCESS) {
        free(coeffs);
        free(paulis);
        return ret_val;
    }

    gh->coeffs = coeffs;
    gh->paulis = paulis;

    return EXIT_SUCCESS;
}


void
h5_grp_hamiltonian_free(h5_grp_hamiltonian gh) {
    free(gh.coeffs);
    free(gh.paulis);
}