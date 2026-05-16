#ifndef CIRC_H
#define CIRC_H

#include <stddef.h>
#include <stdint.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/state_prep_coeff.h"

struct circ_hamil {
	uint32_t qb;
	struct circ_hamil_term {
		double cf;
		struct paulis op;
	} *terms;
	size_t len;
};

struct circ_muldet {
	struct {
		uint64_t idx;
		_Complex double cf;
	} *dets;
	size_t len;
};

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

struct circ {
	struct circ_hamil hm;
	struct circ_muldet md;
	struct circ_cache *cache;
	struct circ_values vals;
	struct qreg reg;
	enum stprep_kind stprep_kind;
	struct data_coeff_matrix cm;
};

int circ_hamil_init(struct circ_hamil *hm, uint32_t qb, size_t len);
void circ_hamil_free(struct circ_hamil *hm);
void circ_hamil_sort_lex(struct circ_hamil *hm);

int circ_muldet_init(struct circ_muldet *md, size_t len);
void circ_muldet_free(struct circ_muldet *md);

/* Release a coefficient-matrix carrier populated by
 * data_coeff_matrix_load() (see ph2run/data.h).  Frees the
 * top-level and per-block C arrays and zeros the struct.
 * Lives here because it is pure pointer cleanup and does
 * not depend on HDF5. */
void data_coeff_matrix_free(struct data_coeff_matrix *cm);

int circ_values_init(struct circ_values *vals, size_t len);
void circ_values_free(struct circ_values *vals);

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
void circ_free(struct circ *ct);
int circ_prepst(struct circ *ct);
int circ_step(struct circ *ct, const struct circ_hamil *hm, double omega);
int circ_step_reverse(
	struct circ *ct, const struct circ_hamil *hm, double omega);
_Complex double circ_measure(struct circ *ct);

#endif // CIRC_H
