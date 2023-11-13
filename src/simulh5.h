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
} simulh5_res;

typedef struct {
    size_t num_qubits;
    size_t num_sum_terms;
    double *coeffs;
    unsigned char *paulis;
} simulh5_grp_pauli_hamil;

simulh5_res
simulh5_grp_pauli_hamil_read(hid_t group_id, simulh5_grp_pauli_hamil *);

void
simulh5_grp_pauli_hamil_drop(simulh5_grp_pauli_hamil);

typedef struct {
    size_t num_steps;
    double *times;
    double *values_real;
    double *values_imag;
} simulh5_grp_time_series;

simulh5_res
simulh5_grp_time_series_read_times(hid_t group_id, simulh5_grp_time_series *);

simulh5_res
simulh5_grp_time_series_write_values(hid_t group_id, simulh5_grp_time_series);

void
simulh5_grp_time_series_drop(simulh5_grp_time_series);

#endif //PHASE2_SIMULH5_H
