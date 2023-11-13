#ifndef PHASE2_SIMULH5_H
#define PHASE2_SIMULH5_H

#include <stdlib.h>
#include "hdf5.h"

#define SIMULH5_GRP_PAULI_HAMIL "pauli_hamil"
#define SIMULH5_GRP_PAULI_HAMIL_COEFFS "coeffs"
#define SIMULH5_GRP_PAULI_HAMIL_PAULIS "paulis"

typedef enum {
    SIMULH5_OK,
    SIMULH5_ERR
} simulh5_res;

typedef struct h5_grp_pauli_hamil {
    size_t num_qubits;
    size_t num_sum_terms;
    double *coeffs;
    unsigned char *paulis;
} simulh5_grp_pauli_hamil;

simulh5_res
simulh5_grp_hamiltonian_read(hid_t group_id, simulh5_grp_pauli_hamil *);

void
simulh5_grp_hamiltonian_drop(simulh5_grp_pauli_hamil);

#endif //PHASE2_SIMULH5_H
