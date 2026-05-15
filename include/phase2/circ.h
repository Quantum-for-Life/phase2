#ifndef CIRC_H
#define CIRC_H

#include <stddef.h>
#include <stdint.h>
#include <time.h>

#include "ph2run/data.h"
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

struct circ_prog {
	unsigned pc;	/* percent of len */
	size_t i, len;
	struct timespec t0;	/* wall-clock start for ETA */
	const char *unit;	/* "step", "sample"; never NULL */
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

void circ_prog_init(struct circ_prog *prog, size_t len, const char *unit);
void circ_prog_tick(struct circ_prog *prog);
/* Format a rich progress line and emit it at info level
 * through the given subsystem tag:
 *
 *   step 17/100 (17%) elapsed 4.31s eta 21.0s
 *
 * Callers control cadence; the library does not throttle. */
void circ_prog_emit(const struct circ_prog *prog, const char *subsys);

int circ_values_init(struct circ_values *vals, size_t len);
void circ_values_free(struct circ_values *vals);

int circ_init(struct circ *ct, data_id fid, size_t vals_len);
void circ_free(struct circ *ct);
int circ_prepst(struct circ *ct);
int circ_step(struct circ *ct, const struct circ_hamil *hm, double omega);
int circ_step_reverse(
	struct circ *ct, const struct circ_hamil *hm, double omega);
_Complex double circ_measure(struct circ *ct);

#endif // CIRC_H
