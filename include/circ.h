/**
 * Circuit interface.
 */
#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

#include "data2.h"

typedef size_t qbid;

struct circ;

struct circuit {
	const char *name;

	size_t num_mea_qb;
	size_t num_sys_qb;
	size_t num_anc_qb;

	int (*reset)(struct circ *);

	int (*prepst)(struct circ *);
	int (*effect)(struct circ *);
	int (*measure)(struct circ *);
};

/*
 * Initialize global circuit environment.
 *
 * Although this function is exposed in the public API, the user doesn't need
 * to call it directly.  It will be called by circ_create() as well.
 *
 * Similarly, the function circ_shutdown() belows will be scheduled with
 * atexit() to clear the global environment during program exit.
 *
 * The initialization may take some time, especially when MPI mode is enabled.
 * It the user cannot wait longer than usual during the first call
 * to circ_init(), it is recomended that that circ_initialize() be called
 * directly at the beginning of the program.
 *
 * This function can be called multiple times from different threads.
 *
 * Return value:	 0	if the environment was successfully initialized
 *				by this function call
 *			 1	if the environment has already been initialized
 *				by another call to this function
 *			-1	in case of failure
 */
int circ_initialize(void);

/*
 * Shut down global environment.
 *
 * This function is registered (using atexit()) during the first call to
 * circ_initialize() to be called automatically when the application closes. The
 * user doesn't need to call this function directly.
 */
void circ_shutdown(void);

/*
 * Create instance of a circuit.
 *
 * Using the circuit description `ct`, allocate memory for a new circuit
 * instance.  The pointer `data` will be stored and will be available during
 * call to `circ_simulate()`.
 */
struct circ *circ_create(struct circuit *ct, void *data);

/*
 * Destroy circuit instance.
 *
 * Free allocated memory.  After the function call, the pointer `c` will no
 * longer point to a valid structure and must be discarded.
 */
void circ_destroy(struct circ *c);

/*
 * Retrieve generic pointer to user data.
 */
void *circ_data(const struct circ *c);

/*
 * Print information about the circuit instance to standard output.
 *
 * Return value:	 0	if successful
 *			-1	in case of failure
 */
int circ_report(struct circ const *c);

/*
 * Reset circuit.
 *
 * Set classical register to zero and call `reset()` function specified by
 * the circuit template (if any).
 *
 * Return value:	 0	if successful
 *			-1	in case of failure
 */
int circ_reset(struct circ *c);

/*
 * Run circuit simulation.
 *
 * This will reset the circuit with circ_reset() and call functions specified by
 * the circuit template in the following order:
 *	prepst(),
 *	effect(),
 *	measure()
 *
 * Return value:	 0	if successful
 *			-1	in case of failure
 */
int circ_run(struct circ *c);

/*
 * Return number of qubits in the "measurement" register.
 */
size_t circ_num_meaqb(const struct circ *c);

/*
 * Return number of qubits in the "system" register.
 */
size_t circ_num_sysqb(const struct circ *c);

/*
 * Return number of qubits in the "ancilla" register.
 */
size_t circ_num_ancqb(const struct circ *c);

qbid circ_meaqb(const struct circ *c, size_t idx);

qbid circ_sysqb(const struct circ *c, size_t idx);

qbid circ_ancqb(const struct circ *c, size_t idx);

void circ_ops_hadamard(struct circ *c, qbid qb);

void circ_ops_sgate(struct circ *c, qbid qb);

double circ_ops_prob0(struct circ *c, qbid qb);

void circ_ops_blankstate(struct circ *c);

void circ_ops_setsysamp(struct circ *c, size_t idx, _Complex double amp);

void circ_ops_crotpauli(struct circ *c, int *paulis, double angle);

typedef unsigned long pauli_pak_t;

/*
 * Hamiltonian module.
 */

/*
 * Structure representing a Hamiltonian as a real linear combination of Pauli
 * strings (codes).  The Pauli operators are packed as a continuous string
 * of pairs of bits.
 */
struct circ_hamil {
	size_t	     num_qubits; /* number of qubits */
	size_t	     num_terms; /* number of terms in the sum */
	double      *coeffs; /* array of coefficients */
	pauli_pak_t *pak; /* array of Pauli operators */
};

/*
 * Initialize Hamiltonian.
 *
 * This function must be called first before any other operation on the
 * structure is performed.
 */
void circ_hamil_init(struct circ_hamil *h);

/*
 * Destroy Hamiltonian.
 *
 * Free allocated memory.
 */
void circ_hamil_destroy(struct circ_hamil *h);

/*
 * Parse data from the open file represented by `fid` descriptor.
 *
 * Return value:  	 0	Data was parsed successfully
 *			-1	Error while reading data
 */
int circ_hamil_from_data2(struct circ_hamil *h, data2_id fid);

/*
 * Retrieve a single Pauli string corresponding to the index `n` in the sum of
 * terms.
 *
 * The value of `n` must be smaller than `h->num_qubits`.  The array `paulis`
 * must be of size at least `n`.
 *
 * After the call to this function, the value of `paulis` will be overwritten
 * with numbers specifying operators in the Pauli string:
 *	0	- I
 *	1	- X
 *	2	- Y
 *	3	- Z
 */
void circ_hamil_paulistr(const struct circ_hamil *h, size_t n, int *paulis);

#endif // PHASE2_CIRC_H
