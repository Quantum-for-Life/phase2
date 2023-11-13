//
// Created by mm on 11/12/23.
//

#include "hdf5.h"

#include "simulh5.h"
#include "log/log.h"

simulh5_res
simulh5_grp_hamiltonian_read(hid_t group_id, simulh5_grp_pauli_hamil *ph) {

    log_debug("Read group " SIMULH5_GRP_PAULI_HAMIL);
    hid_t dset_coeffs_id = H5Dopen2(group_id,
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

    hid_t dset_paulis_id = H5Dopen2(group_id, SIMULH5_GRP_PAULI_HAMIL_PAULIS,
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

void simulh5_grp_hamiltonian_drop(simulh5_grp_pauli_hamil ph) {
    free(ph.coeffs);
    free(ph.paulis);
}
