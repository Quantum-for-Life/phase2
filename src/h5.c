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
    hid_t attr;
    int ret_val = EXIT_SUCCESS;

    if ((attr = H5Aopen(obj_id, attr_name, H5P_DEFAULT)) == H5I_INVALID_HID) {
        ret_val = EXIT_FAILURE;
        goto fail_attr;
    }
    // read the attribute value
    if (H5Aread(attr, type_id, value) < 0) {
        ret_val = EXIT_FAILURE;
    }

    H5Aclose(attr);
    fail_attr:
    return ret_val;
}


h5_group_hamiltonian *
h5_parse_group_hamiltonian(hid_t group_id) {
    h5_group_hamiltonian gh;

    h5_read_attr(group_id, "num_qubits", H5T_NATIVE_ULONG, &gh.num_qubits);
    h5_read_attr(group_id, "num_sum_terms", H5T_NATIVE_ULONG,
                 &gh.num_sum_terms);

    double *coeffs = malloc(sizeof(double) * gh.num_qubits);
    if (coeffs == NULL) {
        return NULL;
    }

    return NULL;
}