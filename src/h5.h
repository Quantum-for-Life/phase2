#ifndef PHASE2_H5_H
#define PHASE2_H5_H

#include <stdlib.h>
#include "hdf5.h"

typedef struct h5_group_hamiltonian_ {
    size_t num_qubits;
    size_t num_sum_terms;
    double *coeffs;
    unsigned char *paulis;
} h5_group_hamiltonian;

typedef struct h5_group_time_series_ {
    size_t steps;
    double *times;
    double *values_real;
    double *values_imag;
} h5_group_time_series;

h5_group_hamiltonian *
h5_parse_group_hamiltonian(hid_t group_id);

#endif //PHASE2_H5_H
