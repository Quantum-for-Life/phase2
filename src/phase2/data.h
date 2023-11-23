#ifndef PHASE2_DATA_H
#define PHASE2_DATA_H

#include <stdlib.h>
#include "hdf5.h"

enum
{
        DATA_OK = 0,
        DATA_ERR,
};

#define DATA_PAULI_HAMIL "pauli_hamil"
#define DATA_PAULI_HAMIL_COEFFS "coeffs"
#define DATA_PAULI_HAMIL_PAULIS "paulis"

struct data_pauli_hamil
{
        size_t num_qubits;
        size_t num_sum_terms;
        double* coeffs;
        unsigned char* paulis;
};

void data_pauli_hamil_init(struct data_pauli_hamil*);

void data_pauli_hamil_destroy(struct data_pauli_hamil*);

int data_pauli_hamil_read(struct data_pauli_hamil*, hid_t);


#define DATA_TIME_SERIES "time_series"
#define DATA_TIME_SERIES_TIMES "times"
#define DATA_TIME_SERIES_VALUES "values"

struct data_time_series
{
        size_t num_steps;
        double* times;
        double* values;
};

void data_time_series_init(struct data_time_series*);

void data_time_series_destroy(struct data_time_series*);

int data_time_series_read(struct data_time_series*, hid_t);

int data_time_series_write(struct data_time_series, hid_t);

#endif //PHASE2_DATA_H
