#ifndef CIRC_H
#define CIRC_H

/*
 * phase2/circ.h -- circuit context, data carriers, and
 * the lifecycle / step / measure API.
 *
 * `struct circ` is the per-simulation aggregate: it owns
 * a Pauli Hamiltonian, a state-prep payload (multidet or
 * coeff_matrix), an MPI-distributed register, a batch
 * cache, and the per-step results buffer.
 *
 * The carriers below (`circ_hamil`, `circ_muldet`,
 * `circ_values`, plus the coefficient-matrix payloads
 * `data_coeff_matrix` / `data_coeff_block`) are populated
 * either by the data subsystem (`ph2run/data.h`) or in-
 * tree by test code.  The free helpers here are pure
 * pointer cleanup and do not depend on HDF5, so circ /
 * libphase2.so can release them without linking the I/O
 * layer.
 *
 * Algorithms (`circ/trott.c`, `trott2.c`, `qdrift.c`,
 * `cmpsit.c`) build a `struct circ`, prepare the state,
 * run a sequence of `circ_step` / `circ_step_reverse`
 * calls, and read out each step's overlap with
 * `circ_measure`.
 */

#include <stddef.h>
#include <stdint.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/state_prep_coeff.h"

/*
 * Sparse Hamiltonian: a sum of weighted Pauli strings.
 *
 *   qb     qubit count the operators apply to.
 *   terms  array of (coefficient, packed Pauli) pairs.
 *   len    number of terms.
 *
 * Allocated by `circ_hamil_init`; freed by
 * `circ_hamil_free`.  Ownership is transferred to a
 * `struct circ` on a successful `circ_init`.
 */
struct circ_hamil {
	uint32_t qb;
	struct circ_hamil_term {
		double cf;
		struct paulis op;
	} *terms;
	size_t len;
};

/*
 * Multi-determinant state-prep payload: a sparse
 * superposition |psi> = sum_i cf_i |idx_i>.
 *
 *   dets   array of (basis-state index, complex
 *          amplitude) pairs.
 *   len    number of basis states.
 */
struct circ_muldet {
	struct {
		uint64_t idx;
		_Complex double cf;
	} *dets;
	size_t len;
};

/*
 * Per-step results: one complex overlap per Trotter
 * step (or qDRIFT / composite sample), filled in by
 * the algorithm during simulation.
 */
struct circ_values {
	_Complex double *z;
	size_t len;
};

/*
 * State-prep subtypes carried in simul.h5.
 *
 * Exactly one of /state_prep/multidet or /state_prep/coeff_matrix
 * must be present; both-present is an error, neither-present is
 * an error.  See data_state_prep_kind() in ph2run/data.h.
 */
enum stprep_kind {
	STPREP_MULTIDET = 1,
	STPREP_COEFF_MATRIX = 2,
};

/*
 * Raw coefficient-matrix data read from
 * /state_prep/coeff_matrix.
 *
 *   nqb           total qubit count (== n_qubits attribute)
 *   n_sites       spatial-orbital count
 *   n_alpha       alpha-spin occupation
 *   n_beta        beta-spin occupation
 *   closed_shell  0 or 1; 1 => C_beta buffers absent
 *   tapered       0 or 1; 1 => bits 0 and n_sites are dropped
 *                              per generated bitstring
 *
 *   n_components  CSF superposition arity.  0 => the trial
 *                 state is a single block carried by the
 *                 top-level C_alpha / C_beta arrays;
 *                 blocks == NULL.  > 0 => the top-level
 *                 C_alpha / C_beta are NULL and `blocks[]`
 *                 carries n_components entries, each with
 *                 its own weight `cf` and C_alpha / C_beta.
 *   C_alpha       single-block alpha coefficients,
 *                 n_sites * n_alpha doubles row-major, or
 *                 NULL when n_components > 0.
 *   C_beta        single-block beta coefficients,
 *                 n_sites * n_beta doubles row-major, or
 *                 NULL when closed_shell or n_components > 0.
 *   blocks        CSF blocks, n_components entries, or NULL
 *                 for the single-block case.  Each block's
 *                 C_alpha / C_beta have the same shapes as
 *                 the top-level arrays.
 *
 * Loaded by data_coeff_matrix_load() (see ph2run/data.h) and
 * released by data_coeff_matrix_free().
 */
struct data_coeff_block {
	double cf;
	const double *C_alpha;
	const double *C_beta;
};

struct data_coeff_matrix {
	uint32_t nqb;
	uint32_t n_sites;
	uint32_t n_alpha;
	uint32_t n_beta;
	int closed_shell;
	int tapered;
	size_t n_components;
	const double *C_alpha;
	const double *C_beta;
	struct data_coeff_block *blocks;
};

/*
 * Per-simulation aggregate.  Members:
 *
 *   hm           sparse Pauli Hamiltonian; adopted by
 *                circ_init, freed by circ_free.
 *   md           multidet state-prep payload; populated
 *                when stprep_kind == STPREP_MULTIDET.
 *   cache        Pauli-rotation batch cache (private
 *                type, allocated in circ_init, freed in
 *                circ_free).
 *   vals         per-step results buffer.
 *   reg          MPI-distributed quantum register.
 *   stprep_kind  selects which state-prep payload to
 *                consume (md vs cm).
 *   cm           coefficient-matrix state-prep payload;
 *                populated when stprep_kind ==
 *                STPREP_COEFF_MATRIX.
 *
 * Callers do not poke at internal fields directly --
 * the lifecycle goes through circ_init / circ_free,
 * the simulation through circ_prepst / circ_step /
 * circ_measure.
 */
struct circ {
	struct circ_hamil hm;
	struct circ_muldet md;
	struct circ_cache *cache;
	struct circ_values vals;
	struct qreg reg;
	enum stprep_kind stprep_kind;
	struct data_coeff_matrix cm;
	struct state_prep_coeff_scratch sp_scratch;
};


/* -- data carriers ------------------------------------------------------- */

/*
 * circ_hamil_init / _free / _sort_lex -- allocate,
 * release, and pre-sort a Hamiltonian carrier.
 *
 * `_init` allocates `len` term slots and records the
 * qubit count; returns 0 on success, -1 on
 * allocation failure.
 *
 * `_free` releases the term buffer and zeroes both the
 * pointer and `len`, so a second free is a clean
 * no-op.  Safe to call on a default-zero `struct
 * circ_hamil`.
 *
 * `_sort_lex` reorders terms by lexicographic Pauli
 * code, which groups terms sharing an MPI exchange
 * contiguously and maximises the batch cache hit rate
 * in circ_step.  Cheap to run; do once before the
 * first step.
 */
int circ_hamil_init(struct circ_hamil *hm, uint32_t qb, size_t len);
void circ_hamil_free(struct circ_hamil *hm);
void circ_hamil_sort_lex(struct circ_hamil *hm);

/*
 * circ_muldet_init / _free -- allocate and release the
 * multi-determinant state-prep carrier.  `_init`
 * returns 0 on success, -1 on allocation failure;
 * `_free` is idempotent.
 */
int circ_muldet_init(struct circ_muldet *md, size_t len);
void circ_muldet_free(struct circ_muldet *md);

/* Release a coefficient-matrix carrier populated by
 * data_coeff_matrix_load() (see ph2run/data.h).  Frees the
 * top-level and per-block C arrays and zeros the struct.
 * Lives here because it is pure pointer cleanup and does
 * not depend on HDF5. */
void data_coeff_matrix_free(struct data_coeff_matrix *cm);

/*
 * circ_values_init / _free -- allocate and release the
 * per-step results buffer (`len` complex slots).
 * `_init` returns 0 on success, -1 on allocation
 * failure; `_free` is idempotent.
 */
int circ_values_init(struct circ_values *vals, size_t len);
void circ_values_free(struct circ_values *vals);


/* -- lifecycle ----------------------------------------------------------- */

/*
 * circ_init - adopt pre-loaded Hamiltonian and state-prep
 * data into a fresh circ context.
 *
 *   hm        Pauli-Hamiltonian, packed (cf, struct paulis).
 *             Passed by value; circ takes ownership of the
 *             buffers and frees them in circ_free.  The
 *             caller must not free them after a successful
 *             circ_init.
 *   sp_kind   selects which state-prep payload to adopt.
 *   sp_data   pointer to a struct circ_muldet (when sp_kind
 *             == STPREP_MULTIDET) or struct data_coeff_matrix
 *             (when STPREP_COEFF_MATRIX).  Copied by value;
 *             same ownership transfer as hm.
 *   vals_len  number of per-step result slots to allocate.
 *
 * Returns 0 on success, -1 on error (with log_error).  On
 * error the function frees any data it has adopted; the
 * caller's locals must not be freed a second time.
 */
int circ_init(struct circ *ct, struct circ_hamil hm,
	enum stprep_kind sp_kind, const void *sp_data, size_t vals_len);

/* Release every owned buffer (Hamiltonian, state-prep
 * payload, cache, register, results) and zero the
 * cache pointer.  Safe on a fully-initialised circ. */
void circ_free(struct circ *ct);


/* -- simulation surface -------------------------------------------------- */

/*
 * circ_prepst -- write the initial state encoded by the
 * adopted state-prep payload into the MPI register.
 * Zeroes the register first.  Returns 0 on success, -1
 * on error.
 */
int circ_prepst(struct circ *ct);

/*
 * circ_step / circ_step_reverse -- apply one Trotter step
 * exp(i*omega*H), traversing the Hamiltonian forward or
 * backward respectively.  Forward and reverse sweeps with
 * the same omega cancel only when each Pauli rotation
 * commutes with the next, hence the order matters for
 * Strang-style 2nd-order Trotter.
 *
 * Returns 0 on success, -1 on error (e.g. a cache flush
 * could not place the next term).
 */
int circ_step(struct circ *ct, const struct circ_hamil *hm, double omega);
int circ_step_reverse(
	struct circ *ct, const struct circ_hamil *hm, double omega);

/*
 * circ_measure -- compute <trial | evolved>, where the
 * trial state is determined by the adopted state-prep
 * payload and `evolved` is the register's current
 * contents.  Returns the complex overlap.
 *
 * Sentinel: returns 0+0i if `ct->stprep_kind` is not
 * one of the known enum values (the switch has no
 * default; this path is unreachable on a well-formed
 * circ but kept defensively).
 */
_Complex double circ_measure(struct circ *ct);

#endif // CIRC_H
