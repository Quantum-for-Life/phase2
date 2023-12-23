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
 * Arguments:
 *	filename	Path to the data file.
 *
 * Return value:
 *	A valid data_id value or DATA_INVALID_FID in case of error
 */
data_id
data2_open(const char *filename);

/**
 * Close data file.
 *
 * Close the file that has been open with the call to data2_open().  Free the
 * resources.
 */
void data2_close(data_id);

/**
 * Get the number of qubits and terms for the "multidet" group.
 *
 * After a successful call, the value pointed to by the argument "num_qubits"
 * stores the number of qubits saved for the "multidet" group.  Similarly, the
 * argument "num_dets" stores the number of terms in the multidet
 * representation.  If the value cannot be read, the function returns '-1',
 * an the value variables poited to is unchanged.
 *
 * Arguments:
 *  fid			Open file id obtained from data2_open()
 *  num_qubits
 *  num_terms	Pointer to a variable where the result will be stored.
 *
 *
 * Return value:
 *   0			if the value was successfully retrieved
 *  -1			in case of error
 */
int
data2_multidet_getnums(data_id fid, size_t *num_qubits, size_t *num_dets);

/**
 * Perform action "op" on each determinant in multidet group.
 *
 * Call user-supplied "op" function on each Slater determinant in the
 * "multidet" group.  The operator "op" takes as an argument a complex
 * coeffitient for the determinant and the "index" of the determinant
 * (given the qubit values in the binary reresentation, least significant bit
 * represented by qubit 0); and a generic pointer to the data specified
 * by the user.
 *
 * The return value of the operator "op" controls the iteration.  If "op"
 * returns "0", the iteration will continue with the next element.  If the
 *return value is non-zero, the iteration will stop and the value is returned to
 *the caller.  By convention, a negative value should indicate an error, whereas
 *a positive value should mean that the iteration terminated early but with no
 *error.
 *
 * Return value:
 *   0      if the full iteration completed sucessfully
 *  -1      if the data could not be retrieved,
 *  or a user-defined value, if the iteration was terminated early
 */
int
data2_multidet_foreach(
	data_id fid, int (*op)(_Complex double, size_t, void *), void *op_data);

int
data2_hamil_getnums(data_id fid, size_t *num_qubits, size_t *num_terms);

int
data2_hamil_foreach(data_id fid,
	int (*op)(_Complex double, unsigned char *, void *), void *op_data);

int
data2_times_getnums(data_id fid, size_t &num_steps);

int
data2_times_foreach(
	data_id fid, int (*op)(double, _Complex double, void *), void *op_data);


/* ---------------------------------------------------------------------------
 * This is a deprecated API.  To be removed.
 */

struct data_state_prep_multidet {
	size_t		 num_qubits;
	size_t		 num_terms;
	_Complex double *coeffs;
	unsigned char	*dets;
};

struct data_state_prep {
	struct data_state_prep_multidet multidet;
};

struct data_pauli_hamil {
	size_t	       num_qubits;
	size_t	       num_terms;
	double	      *coeffs;
	unsigned char *paulis;
	double	       norm;
};

struct data_time_series {
	size_t		 num_steps;
	double		*times;
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

/* ---------------------------------------------------------------------------*/

#endif // DATA_H
