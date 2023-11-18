#ifndef PHASE2_SDAT_H
#define PHASE2_SDAT_H

#include <stdlib.h>
#include "hdf5.h"

enum {
    sdat_ok,
    sdat_err
};

#define SDAT_PAULI_HAMIL "pauli_hamil"
#define SDAT_PAULI_HAMIL_COEFFS "coeffs"
#define SDAT_PAULI_HAMIL_PAULIS "paulis"

typedef struct {
    size_t num_qubits;
    size_t num_sum_terms;
    double *coeffs;
    unsigned char *paulis;
} sdat_pauli_hamil;

void
sdat_pauli_hamil_init(sdat_pauli_hamil *);

void
sdat_pauli_hamil_drop(sdat_pauli_hamil);

int
sdat_pauli_hamil_read(sdat_pauli_hamil *, hid_t);


#define SDAT_TIME_SERIES "time_series"
#define SDAT_TIME_SERIES_TIMES "times"
#define SDAT_TIME_SERIES_VALUES "values"

typedef struct {
    size_t num_steps;
    double *times;
    double *values;
} sdat_time_series;

void
sdat_time_series_init(sdat_time_series *);

void
sdat_time_series_drop(sdat_time_series);

int
sdat_time_series_read(sdat_time_series *, hid_t);

int
sdat_time_series_write(sdat_time_series, hid_t);

#endif //PHASE2_SDAT_H
