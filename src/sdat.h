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

struct sdat_pauli_hamil {
    size_t num_qubits;
    size_t num_sum_terms;
    double *coeffs;
    unsigned char *paulis;
};

void
sdat_pauli_hamil_init(struct sdat_pauli_hamil *);

void
sdat_pauli_hamil_destroy(struct sdat_pauli_hamil *);

int
sdat_pauli_hamil_read(struct sdat_pauli_hamil *, hid_t);


#define SDAT_TIME_SERIES "time_series"
#define SDAT_TIME_SERIES_TIMES "times"
#define SDAT_TIME_SERIES_VALUES "values"

struct sdat_time_series {
    size_t num_steps;
    double *times;
    double *values;
};

void
sdat_time_series_init(struct sdat_time_series *);

void
sdat_time_series_destroy(struct sdat_time_series *);

int
sdat_time_series_read(struct sdat_time_series *, hid_t);

int
sdat_time_series_write(struct sdat_time_series, hid_t);

#endif //PHASE2_SDAT_H
