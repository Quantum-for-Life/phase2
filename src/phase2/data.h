#ifndef PHASE2_DATA_H
#define PHASE2_DATA_H

#include <stdlib.h>

#define DATA_INVALID_OBJID  (-1)

#define DATA_STATE_PREP "state_prep"

#define DATA_STATE_PREP_MULTIDET "multidet"
#define DATA_STATE_PREP_MULTIDET_COEFFS "coeffs"
#define DATA_STATE_PREP_MULTIDET_DETS "dets"

#define DATA_PAULI_HAMIL "pauli_hamil"
#define DATA_PAULI_HAMIL_COEFFS "coeffs"
#define DATA_PAULI_HAMIL_PAULIS "paulis"
#define DATA_PAULI_HAMIL_NORM "normalization"

#define DATA_TIME_SERIES "time_series"
#define DATA_TIME_SERIES_TIMES "times"
#define DATA_TIME_SERIES_VALUES "values"

enum {
        // Unspecified error
        DATA_ERR = -1,
        // Success
        DATA_OK = 0,
};


typedef int64_t dataid_t; // This is the same as HDF5's hid_t

dataid_t data_file_open(const char *filename);

void data_file_close(dataid_t file_id);

struct data {
        struct data_state_prep *state_prep;
        struct data_pauli_hamil *pauli_hamil;
        struct data_time_series *time_series;
};

void data_init(struct data *dat);

void data_destroy(struct data *dat);

int data_parse(struct data *dat, dataid_t obj_id);


struct data_state_prep {
        struct data_state_prep_multidet *multidet;
};

void data_state_prep_init(struct data_state_prep *dat);

void data_state_prep_destroy(struct data_state_prep *dat);

int data_state_prep_parse(struct data_state_prep *dat, dataid_t obj_id);


struct data_state_prep_multidet {
        size_t num_qubits;
        size_t num_terms;
        double *coeffs;
        unsigned char *dets;
};

void data_state_prep_multidet_init(struct data_state_prep_multidet *dat);

void data_state_prep_multidet_destroy(struct data_state_prep_multidet *dat);

int data_state_prep_multidet_parse(struct data_state_prep_multidet *dat,
                                   dataid_t obj_id);


struct data_pauli_hamil {
        size_t num_qubits;
        size_t num_terms;
        double *coeffs;
        unsigned char *paulis;
        double norm;
};

void data_pauli_hamil_init(struct data_pauli_hamil *dat);

void data_pauli_hamil_destroy(struct data_pauli_hamil *dat);

int data_pauli_hamil_parse(struct data_pauli_hamil *, dataid_t obj_id);


struct data_time_series {
        size_t num_steps;
        double *times;
        double *values;
};

void data_time_series_init(struct data_time_series *dat);

void data_time_series_destroy(struct data_time_series *dat);

int data_time_series_parse(struct data_time_series *dat, dataid_t obj_id);

int data_time_series_write(const struct data_time_series *dat, dataid_t obj_id);

#endif //PHASE2_DATA_H
