#ifndef DATA_H
#define DATA_H

#include <stdlib.h>

#define DATA_INVALID_FID (-1)

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

typedef int64_t data_id; // This is the same as HDF5's hid_t

struct data_state_prep_multidet {
	size_t		 num_qubits;
	size_t		 num_terms;
	_Complex double *coeffs;
	unsigned char   *dets;
};

struct data_state_prep {
	struct data_state_prep_multidet multidet;
};

struct data_pauli_hamil {
	size_t	       num_qubits;
	size_t	       num_terms;
	double	       *coeffs;
	unsigned char *paulis;
	double	       norm;
};

struct data_time_series {
	size_t		 num_steps;
	double	       *times;
	_Complex double *values;
};

struct data {
	struct data_state_prep	state_prep;
	struct data_pauli_hamil pauli_hamil;
	struct data_time_series time_series;
};

data_id
data_file_open(const char *filename);

void
data_file_close(data_id fid);

void
data_init(struct data *dat);

void
data_destroy(struct data *dat);

int
data_parse(struct data *dat, data_id fid);

int
data_time_series_write(data_id fid, const struct data_time_series *dat);

#endif // DATA_H
