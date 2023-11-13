#ifndef PHASE2_H5_H
#define PHASE2_H5_H

#include <stdlib.h>
#include "hdf5.h"

#define H5_GRP_HAMILTONIAN "hamiltonian"

typedef struct h5_grp_hamiltonian {
    int num_qubits;
    int num_sum_terms;
    double *coeffs;
    int *paulis;
} h5_grp_hamiltonian;

int
h5_grp_hamiltonian_read(hid_t group_id, h5_grp_hamiltonian *);

void
h5_grp_hamiltonian_free(h5_grp_hamiltonian);

#endif //PHASE2_H5_H
