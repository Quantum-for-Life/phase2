#ifndef PHASE2_DATA_H
#define PHASE2_DATA_H

#include <stddef.h>
#include <stdint.h>

#define DATA_INVALID_FID (-1)

/* Group, dataset names */
#define DATA_STPREP "state_prep"
#define DATA_STPREP_MULTIDET "multidet"
#define DATA_STPREP_MULTIDET_COEFFS "coeffs"
#define DATA_STPREP_MULTIDET_DETS "dets"

#define DATA_HAMIL "pauli_hamil"
#define DATA_HAMIL_COEFFS "coeffs"
#define DATA_HAMIL_PAULIS "paulis"
#define DATA_HAMIL_NORM "normalization"

#define DATA_CIRCTROTT "circ_trott"
#define DATA_CIRCTROTT_DELTA "delta"
#define DATA_CIRCTROTT_VALUES "values"

#define DATA_CIRCQDRIFT "circ_qdrift"
#define DATA_CIRCQDRIFT_STEPSIZE "step_size"
#define DATA_CIRCQDRIFT_NUMSAMPLES "num_samples"
#define DATA_CIRCQDRIFT_DEPTH "depth"
#define DATA_CIRCQDRIFT_VALUES "values"

/**
 * Handle to a data file
 */
typedef int64_t data_id;

/**
 * Open data file.
 *
 * Arguments:
 *
 *	filename	Path to the data file.
 *
 * Return value:
 *
 *	A valid data_id value or DATA_INVALID_FID in case of error
 */
data_id data_open(const char *filename);

/**
 * Close data file.
 *
 * Close the file that has been open with the call to data_open().  Free the
 * resources.
 */
void data_close(data_id);

int data_grp_create(data_id fid, const char *grp_name);

#define data_attr_read(fid, grp_name, attr_name, attr_buf)                     \
	_Generic((attr_buf), double *: data_attr_read_dbl)(                    \
		fid, grp_name, attr_name, attr_buf)

#define data_attr_write(fid, grp_name, attr_name, attr)                        \
	_Generic((attr), double: data_attr_write_dbl)(                         \
		fid, grp_name, attr_name, attr)

int data_attr_read_dbl(
	data_id fid, const char *grp_name, const char *attr_name, double *a);

/* Create an attribute and write the value of 'a' to it. */
int data_attr_write_dbl(
	data_id fid, const char *grp_name, const char *attr_name, double a);

/**
 * Get the number of qubits and terms for the "multidet" group.
 *
 * After a successful call, the value pointed to by the argument "num_qubits"
 * stores the number of qubits saved for the "multidet" group.  Similarly, the
 * argument "num_dets" stores the number of terms in the multidet
 * representation.  If the value cannot be read, the function returns '-1', and
 * the value of the variables poited to remains unchanged.
 *
 * Arguments:
 *
 *  fid			Open file id obtained from data_open()
 *  num_qubits
 *  num_dets	Pointer to a variable where the result will be stored.
 *
 * Return value:
 *
 *   0			if the value was successfully retrieved
 *  -1			in case of error
 */
int data_multidet_getnums(data_id fid, uint32_t *nqb, size_t *ndets);

/**
 * Perform action "op" on each determinant in "multidet" group.
 *
 * Call user-supplied "op" function on each Slater determinant in the
 * "multidet" group.  The operator "op" takes as arguments a complex
 * coeffitient and the "index" of the determinant (given the qubit values in the
 * binary reresentation, the least significant bit represented by qubit 0); and
 * a generic pointer to the data specified by the user.
 *
 * The return value of the operator "op" controls the iteration.  If "op"
 * returns "0", the iteration will continue with the next element.  If the
 * return value is non-zero, the iteration will stop and the value is returned
 * to the caller.  By convention, a negative value should indicate an error,
 * whereas a positive value means that the iteration simply terminated early no
 * error.
 *
 * Return value:
 *
 *   0      if the full iteration completed sucessfully
 *  -1      if the data could not be retrieved,
 *  or a user-defined value, if the iteration was terminated early
 */
int data_multidet_foreach(data_id fid,
	int (*op)(_Complex double cf, uint64_t idx, void *), void *op_data);

/**
 * Get the number of qubits and terms for the "hamil" group.
 *
 * After a successful call, the value pointed to by the argument "num_qubits"
 * stores the number of qubits saved for the "hamil" group.  Similarly, the
 * argument "num_terms" stores the number of terms of the hamiltonian.  If the
 * values cannot be read, the function returns '-1', and the value of the
 * variables poited to remains unchanged.
 *
 * Arguments:
 *
 *  fid			Open file id obtained from data_open()
 *  num_qubits
 *  num_terms	Pointer to a variable where the result will be stored.
 *
 * Return value:
 *
 *   0			if the value was successfully retrieved
 *  -1			in case of error
 */
int data_hamil_getnums(data_id fid, uint32_t *nqb, size_t *nterms);

/**
 * Get the normalization factor for the "hamil" group.
 *
 * After a successful call, the value pointed to by the argument "norm" stores
 * the normalization factor for the "hamil" group.  If the value cannot be read,
 * the function returns '-1', and the value of the variable poited to remains
 * unchanged.
 *
 * Arguments:
 *
 *  fid			Open file id obtained from data_open()
 *  norm		Pointer to a variable where the result will be stored.
 *
 * Return value:
 *
 *   0			if the value was successfully retrieved
 *  -1			in case of error
 */
int data_hamil_getnorm(data_id fid, double *norm);

/**
 * Perform action "op" on each term of the Hamiltonian in "pauli_hamil" group.
 *
 * Call user-supplied "op" function on each term of the Hamiltonian. The
 * operator "op" takes as arguments a real coeffitient and the array of length
 * num_qubits (see function data_hamil_getnums()) filled with values
 * representing Pauli operators:
 *
 *	0 - I (identity)
 *	1 - Pauli X
 *	2 - Pauli Y
 *	3 - Pauli Z
 *
 * and a generic pointer to the data specified by the user.  The value of the
 * array will change with each iteration.  After the iteration finished, the
 * pointer to the array no longer refers to valid memory.  Any attempt to store
 * and use it later is undefined behaviour.
 *
 * The return value of the operator "op" controls the iteration.  If "op"
 * returns "0", the iteration will continue with the next element.  If the
 * return value is non-zero, the iteration will stop and the value is returned
 * to the caller.  By convention, a negative value should indicate an error,
 * whereas a positive value means that the iteration simply terminated early no
 * error.
 *
 * Return value:
 *
 *   0      if the full iteration completed sucessfully
 *  -1      if the data could not be retrieved,
 *  or a user-defined value, if the iteration was terminated early
 */
int data_hamil_foreach(
	data_id fid, int (*op)(double, unsigned char *, void *), void *op_data);

int data_res_write(data_id fid, const char *grp_name, const char *dset_name,
	const _Complex double *vals, size_t nvals);

#endif // PHASE2_DATA_H
