#ifndef PHASE2_RUN_H
#define PHASE2_RUN_H

#include <stddef.h>

/*
 * phase2_run - apply Trotterised Hamiltonian evolution
 * and return the overlap with the initial state.
 *
 * Computes <psi| prod_k exp(i * delta * coeffs[k] * P_k) |psi>
 * where P_k are Pauli strings and |psi> is a computational
 * basis state given by psi_str.
 *
 * The number of qubits is strlen(psi_str).
 *
 * Arguments:
 *   pauli_strs  Array of nterms NUL-terminated Pauli
 *               strings.  Each string is space-separated
 *               tokens like "X0 Z1 Y2".  Unmentioned
 *               qubits are identity.
 *   coeffs      Array of nterms real coefficients.
 *   nterms      Number of Hamiltonian terms.
 *   delta       Scaling factor (time step).
 *   psi_str     Bitstring for the initial state.
 *               Characters must be '0' or '1'.
 *               Length determines nqb.
 *   out_re      Pointer to store the real part.
 *   out_im      Pointer to store the imaginary part.
 *
 * Returns:
 *    0  on success
 *   -1  on error (invalid input, allocation failure)
 *
 * MPI:
 *   Initialises MPI on first call if not already done.
 *   Works without mpirun (single rank).
 */
int phase2_run(const char **pauli_strs,
	const double *coeffs, size_t nterms, double delta,
	const char *psi_str, double *out_re, double *out_im);

#endif /* PHASE2_RUN_H */
