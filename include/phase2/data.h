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

#define DATA_STPREP_COEFFMAT "coeff_matrix"
#define DATA_STPREP_COEFFMAT_CA "C_alpha"
#define DATA_STPREP_COEFFMAT_CB "C_beta"
#define DATA_STPREP_COEFFMAT_NQB "n_qubits"
#define DATA_STPREP_COEFFMAT_NS "n_sites"
#define DATA_STPREP_COEFFMAT_NA "n_alpha"
#define DATA_STPREP_COEFFMAT_NB "n_beta"
#define DATA_STPREP_COEFFMAT_CS "closed_shell"
#define DATA_STPREP_COEFFMAT_TAP "tapered"
#define DATA_STPREP_COEFFMAT_CSF "csf"
#define DATA_STPREP_COEFFMAT_CSF_NCOMP "n_components"
#define DATA_STPREP_COEFFMAT_CSF_CF "coefficient"

#define DATA_HAMIL "pauli_hamil"
#define DATA_HAMIL_COEFFS "coeffs"
#define DATA_HAMIL_NORM "normalization"
#define DATA_HAMIL_PAULIS "paulis"

#define DATA_CIRCTROTT "circ_trott"
#define DATA_CIRCTROTT_DELTA "delta"
#define DATA_CIRCTROTT_VALUES "values"

#define DATA_CIRCTROTT2 "circ_trott2"
#define DATA_CIRCTROTT2_DELTA "delta"
#define DATA_CIRCTROTT2_VALUES "values"

#define DATA_CIRCQDRIFT "circ_qdrift"
#define DATA_CIRCQDRIFT_DEPTH "depth"
#define DATA_CIRCQDRIFT_NUMSAMPLES "num_samples"
#define DATA_CIRCQDRIFT_STEPSIZE "step_size"
#define DATA_CIRCQDRIFT_VALUES "values"
#define DATA_CIRCQDRIFT_SEED "seed"

#define DATA_CIRCCMPSIT "circ_cmpsit"
#define DATA_CIRCCMPSIT_DEPTH "depth"
#define DATA_CIRCCMPSIT_LENGTH "length"
#define DATA_CIRCCMPSIT_ANGLEDET "angle_det"
#define DATA_CIRCCMPSIT_ANGLERAND "angle_rand"
#define DATA_CIRCCMPSIT_STEPS "steps"
#define DATA_CIRCCMPSIT_VALUES "values"
#define DATA_CIRCCMPSIT_SEED "seed"

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

#define DECL_DATA_ATTR_READ(suff, type)                                        \
	int data_attr_read_##suff(data_id fid, const char *grp_name,           \
		const char *attr_name, type *a);

DECL_DATA_ATTR_READ(i, int);
DECL_DATA_ATTR_READ(ul, unsigned long);
DECL_DATA_ATTR_READ(dbl, double);

#define data_attr_read(fid, grp_name, attr_name, attr_buf)                     \
	_Generic((attr_buf),                                                   \
		int *: data_attr_read_i,                                       \
		unsigned long *: data_attr_read_ul,                            \
		double *: data_attr_read_dbl)(                                 \
		fid, grp_name, attr_name, attr_buf)

/* Create an attribute and write the value of 'a' to it. */
#define DECL_DATA_ATTR_WRITE(suff, type)                                       \
	int data_attr_write_##suff(data_id fid, const char *grp_name,          \
		const char *attr_name, type a)

DECL_DATA_ATTR_WRITE(i, int);
DECL_DATA_ATTR_WRITE(ul, unsigned long);
DECL_DATA_ATTR_WRITE(dbl, double);

#define data_attr_write(fid, grp_name, attr_name, attr)                        \
	_Generic((attr),                                                       \
		int: data_attr_write_i,                                        \
		double: data_attr_write_dbl,                                   \
		unsigned long: data_attr_write_ul)(                            \
		fid, grp_name, attr_name, attr)

/**
 * Get the number of qubits and terms for the "multidet" group.
 *
 * After a successful call, the value pointed to by the argument
 * "num_qubits" stores the number of qubits saved for the "multidet"
 * group.  Similarly, the argument "num_dets" stores the number of terms
 * in the multidet representation.  If the value cannot be read, the
 * function returns '-1', and the value of the variables poited to
 * remains unchanged.
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
 * coeffitient and the "index" of the determinant (given the qubit
 * values in the binary reresentation, the least significant bit
 * represented by qubit 0); and a generic pointer to the data specified
 * by the user.
 *
 * The return value of the operator "op" controls the iteration.  If
 * "op" returns "0", the iteration will continue with the next element.
 * If the return value is non-zero, the iteration will stop and the
 * value is returned to the caller.  By convention, a negative value
 * should indicate an error, whereas a positive value means that the
 * iteration simply terminated early no error.
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
 * State-prep subtypes carried in simul.h5.
 *
 * Exactly one of /state_prep/multidet or /state_prep/coeff_matrix
 * must be present; both-present is an error, neither-present is
 * an error.  See data_state_prep_kind().
 */
enum stprep_kind {
	STPREP_MULTIDET = 1,
	STPREP_COEFF_MATRIX = 2,
};

/**
 * Probe simul.h5 for the active state-prep subtype.
 *
 * Inspects the file for /state_prep/multidet and
 * /state_prep/coeff_matrix and applies the dispatch table:
 *
 *  multidet | coeff_matrix | result
 *  ---------+--------------+--------
 *  absent   | absent       | -ENOENT
 *  present  | absent       | STPREP_MULTIDET
 *  absent   | present      | STPREP_COEFF_MATRIX
 *  present  | present      | -EINVAL  (ambiguous; rebuild pak)
 *
 * On success *out holds the selected kind and the function
 * returns 0.  On failure *out is unchanged and the function
 * returns a negative errno-style value.
 *
 * Documented further in phase2/doc/simul-h5-specs.md
 * "dispatch rules".
 */
int data_state_prep_kind(data_id fid, enum stprep_kind *out);

/**
 * Get attributes from /state_prep/coeff_matrix/.
 *
 *   nqb          total qubit count (== n_qubits attribute)
 *   n_sites      spatial-orbital count
 *   n_alpha      alpha-spin occupation
 *   n_beta       beta-spin occupation
 *   closed_shell 0 or 1 (1 => C_beta dataset absent)
 *   tapered      0 or 1 (1 => bits 0 and n_sites are dropped
 *                          per generated bitstring)
 *
 * Returns 0 on success, -1 on error.  Output values are
 * unchanged on error.
 */
int data_coeff_matrix_getnums(data_id fid, uint32_t *nqb, uint32_t *n_sites,
	uint32_t *n_alpha, uint32_t *n_beta, int *closed_shell, int *tapered);

/**
 * Read C_alpha and (optionally) C_beta from the top-level
 * /state_prep/coeff_matrix/ group.
 *
 * Caller must supply buffers sized:
 *   C_alpha: n_sites * n_alpha doubles  (row-major)
 *   C_beta : n_sites * n_beta  doubles  (row-major), or NULL
 *            if closed_shell == 1.
 *
 * Passing a non-NULL C_beta on a closed-shell file or NULL on
 * an open-shell file is a programmer error and returns -1
 * without touching the buffers.
 *
 * Returns 0 on success, -1 on error.
 */
int data_coeff_matrix_read(data_id fid, double *C_alpha, double *C_beta);

/**
 * Count CSF components under /state_prep/coeff_matrix/csf/.
 *
 * Sets *n on success and returns 0.  If the csf/ subgroup is
 * absent, *n is set to 0 and the function returns 0 (the
 * file encodes a single block, not a CSF superposition).
 * Returns -1 on read error (csf/ present but malformed).
 */
int data_coeff_matrix_csf_count(data_id fid, size_t *n);

/**
 * Read CSF component k.
 *
 * Reads /state_prep/coeff_matrix/csf/<k>/coefficient (scalar)
 * and the per-component C_alpha / C_beta datasets.  Buffer
 * shapes match data_coeff_matrix_read().
 *
 * Returns 0 on success, -1 on error.
 */
int data_coeff_matrix_csf_read(data_id fid, size_t k, double *coefficient,
	double *C_alpha, double *C_beta);

/**
 * Get the number of qubits and terms for the "hamil" group.
 *
 * After a successful call, the value pointed to by the argument
 * "num_qubits" stores the number of qubits saved for the "hamil" group.
 * Similarly, the argument "num_terms" stores the number of terms of the
 * hamiltonian.  If the values cannot be read, the function returns
 * '-1', and the value of the variables poited to remains unchanged.
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
 * After a successful call, the value pointed to by the argument "norm"
 * stores the normalization factor for the "hamil" group.  If the value
 * cannot be read, the function returns '-1', and the value of the
 * variable poited to remains unchanged.
 *
 * Arguments:
 *
 *  fid			Open file id obtained from data_open()
 *  norm		Pointer to a variable where the result will be
 * stored.
 *
 * Return value:
 *
 *   0			if the value was successfully retrieved
 *  -1			in case of error
 */
int data_hamil_getnorm(data_id fid, double *norm);

/**
 * Perform action "op" on each term of the Hamiltonian in "pauli_hamil"
 *group.
 *
 * Call user-supplied "op" function on each term of the Hamiltonian. The
 * operator "op" takes as arguments a real coeffitient and the array of
 *length num_qubits (see function data_hamil_getnums()) filled with
 *values representing Pauli operators:
 *
 *	0 - I (identity)
 *	1 - Pauli X
 *	2 - Pauli Y
 *	3 - Pauli Z
 *
 * and a generic pointer to the data specified by the user.  The value
 *of the array will change with each iteration.  After the iteration
 *finished, the pointer to the array no longer refers to valid memory.
 *Any attempt to store and use it later is undefined behaviour.
 *
 * The return value of the operator "op" controls the iteration.  If
 *"op" returns "0", the iteration will continue with the next element.
 *If the return value is non-zero, the iteration will stop and the value
 *is returned to the caller.  By convention, a negative value should
 *indicate an error, whereas a positive value means that the iteration
 *simply terminated early no error.
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
