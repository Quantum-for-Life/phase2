#ifndef DATA_H
#define DATA_H

#include <stdlib.h>

#define DATA_INVALID_FID (-1)

/**
 * Handle to a data file
 */
typedef int64_t data_id;

/**
 * Open data file.
 *
 * Arguments:		filename	Path to the data file.
 *
 * Return value: 	A valid data_id value or DATA_INVALID_FID
 * 			in case of error
 */
data_id
data2_open(const char *filename);

/**
 * Close data file.
 *
 * Close the file that has been open with the call to data2_open().  Free the
 * resources.
 */
void
data2_close(data_id);

/**
 * Get the number of qubits and terms for the "multidet" group.
 *
 * After a successful call, the value pointed to by the argument "num_qubits"
 * stores the number of qubits saved for the "multidet" group.  Similarly, the
 * argument "num_dets" stores the number of terms in the multidet
 * representation.  If the value cannot be read, the function returns '-1',
 * an the value variables poited to is unchanged.
 *
 * Arguments:		fid		Open file id obtained from data2_open()
 * 			num_qubits
 * 			num_terms	Pointer to a variable where the result
 * 					will be stored.
 *
 * Return value:	 0	if the value was successfully retrieved
 * 			-1	in case of error
 */
int
data2_multidet_getnums(data_id fid, size_t *num_qubits, size_t *num_dets);

int
data2_multidet_foreach(
	data_id fid, int (*op)(_Complex double, size_t, void *), void *data);

/* ---------------------------------------------------------------------------
 * This is a deprecated API.  To be removed.
 */

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

struct data_state_prep_multidet {
	size_t           num_qubits;
	size_t           num_terms;
	_Complex double *coeffs;
	unsigned char *  dets;
};

struct data_state_prep {
	struct data_state_prep_multidet multidet;
};

struct data_pauli_hamil {
	size_t         num_qubits;
	size_t         num_terms;
	double *       coeffs;
	unsigned char *paulis;
	double         norm;
};

struct data_time_series {
	size_t           num_steps;
	double *         times;
	_Complex double *values;
};

struct data {
	struct data_state_prep  state_prep;
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

/* ---------------------------------------------------------------------------*/

#endif // DATA_H
