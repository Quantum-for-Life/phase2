#ifndef PHASE2_SIMULH5_H
#define PHASE2_SIMULH5_H

#include <stdlib.h>
#include "hdf5.h"

#define SIMULH5_GRP_PAULI_HAMIL "pauli_hamil"
#define SIMULH5_GRP_PAULI_HAMIL_COEFFS "coeffs"
#define SIMULH5_GRP_PAULI_HAMIL_PAULIS "paulis"

#define SIMULH5_GRP_TIME_SERIES "time_series"
#define SIMULH5_GRP_TIME_SERIES_TIMES "times"
#define SIMULH5_GRP_TIME_SERIES_VALUES_REAL "values_real"
#define SIMULH5_GRP_TIME_SERIES_VALUES_IMAG "values_imag"

typedef enum {
    SIMULH5_OK,
    SIMULH5_ERR
} simulh5_result;

typedef struct {

    struct {
        size_t num_qubits;
        size_t num_sum_terms;
        double *coeffs;
        unsigned char *paulis;
    } pauli_hamil;

    struct {
        size_t num_steps;
        double *times;
        double *values_real;
        double *values_imag;
    } time_series;

} simulh5;

simulh5 *simulh5_create();

simulh5_result
simulh5_read(simulh5 *sh, hid_t obj_id);

void simulh5_free(simulh5 *);

#endif //PHASE2_SIMULH5_H
