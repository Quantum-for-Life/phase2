#ifndef PHASE2_DATA_H
#define PHASE2_DATA_H

#include <stdlib.h>

enum {
        // Unspecified error
        DATA_ERR = -1,
        // Success
        DATA_OK = 0,
};

#define DATA_INVALID_OBJID  (-1)

typedef int64_t dataid_t; // This is the same as HDF5's hid_t

dataid_t data_file_open(const char *filename);

void data_file_close(dataid_t);


#define DATA_PAULI_HAMIL "pauli_hamil"
#define DATA_PAULI_HAMIL_COEFFS "coeffs"
#define DATA_PAULI_HAMIL_PAULIS "paulis"

struct data_pauli_hamil {
        size_t num_qubits;
        size_t num_terms;
        double *coeffs;
        unsigned char *paulis;
};

void data_pauli_hamil_init(struct data_pauli_hamil *);

void data_pauli_hamil_destroy(struct data_pauli_hamil *);

int data_pauli_hamil_read(struct data_pauli_hamil *, dataid_t);


#define DATA_TIME_SERIES "time_series"
#define DATA_TIME_SERIES_TIMES "times"
#define DATA_TIME_SERIES_VALUES "values"

struct data_time_series {
        size_t num_steps;
        double *times;
        double *values;
};

void data_time_series_init(struct data_time_series *);

void data_time_series_destroy(struct data_time_series *);

int data_time_series_read(struct data_time_series *, dataid_t);

int data_time_series_write(const struct data_time_series *, dataid_t);

#endif //PHASE2_DATA_H
