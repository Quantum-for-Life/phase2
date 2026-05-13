#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define LOG_SUBSYS "circ"
#include "log.h"
#include "phase2/circ.h"
#include "phase2/paulis.h"
#include "phase2/state_prep_coeff.h"

#include "circ_cache.h"

int circ_hamil_init(struct circ_hamil *hm, uint32_t qb, size_t len)
{
	hm->terms = malloc(sizeof *hm->terms * len);
	if (!hm->terms)
		return -1;
	hm->len = len;
	hm->qb = qb;

	return 0;
}

void circ_hamil_free(struct circ_hamil *hm)
{
	if (hm->terms != nullptr)
		free(hm->terms);
}

struct hamil_iter_data {
	size_t i;
	double norm;
	struct circ_hamil *hamil;
};

static int hamil_iter(double cf, unsigned char *ops, void *iter_data)
{
	struct hamil_iter_data *id = iter_data;
	struct circ_hamil *h = id->hamil;
	const size_t i = id->i++;
	const uint32_t nqb = h->qb;

	h->terms[i].cf = cf * id->norm;
	struct paulis op = paulis_new();
	for (uint32_t j = 0; j < nqb; j++)
		paulis_set(&op, ops[j], j);
	h->terms[i].op = op;

	return 0;
}

static int circ_hamil_from_file(struct circ_hamil *h, const data_id fid)
{
	uint32_t nqb;
	size_t nterms;
	double norm;

	if (data_hamil_getnums(fid, &nqb, &nterms) < 0)
		return -1;
	if (data_hamil_getnorm(fid, &norm) < 0)
		return -1;
	if (circ_hamil_init(h, nqb, nterms) < 0)
		return -1;

	struct hamil_iter_data id = { .i = 0, .norm = norm, .hamil = h };
	if (data_hamil_foreach(fid, hamil_iter, &id) != 0)
		return -1;

	return 0;
}

static int hamil_term_cmp_lex(const void *a, const void *b)
{
	const struct paulis x = ((const struct circ_hamil_term *)a)->op;
	const struct paulis y = ((const struct circ_hamil_term *)b)->op;

	return paulis_cmp(x, y);
}

/*
 * Lexicographic sort groups terms with the same hi-qubit
 * Pauli code contiguously, maximising cache hits in
 * circ_step_generic and minimising MPI exchanges.
 */
void circ_hamil_sort_lex(struct circ_hamil *hm)
{
	qsort(hm->terms, hm->len, sizeof(struct circ_hamil_term),
		hamil_term_cmp_lex);
	log_debug("hamil_sort_lex: sorted %zu terms", hm->len);
}

int circ_muldet_init(struct circ_muldet *md, size_t len)
{
	md->dets = malloc(sizeof *md->dets * len);
	if (!md->dets)
		return -1;
	md->len = len;

	return 0;
}

void circ_muldet_free(struct circ_muldet *md)
{
	if (md->dets != nullptr)
		free(md->dets);
}

struct iter_muldet_data {
	size_t i;
	struct circ_muldet *muldet;
};

static int iter_muldet(_Complex double cf, const uint64_t idx, void *iter_data)
{
	struct iter_muldet_data *id = iter_data;
	struct circ_muldet *m = id->muldet;
	const size_t i = id->i++;

	m->dets[i].cf = cf;
	m->dets[i].idx = idx;

	return 0;
}

static int circ_muldet_from_file(struct circ_muldet *m, const data_id fid)
{
	uint32_t nqb;
	size_t ndets;

	if (data_multidet_getnums(fid, &nqb, &ndets) < 0)
		return -1;
	if (circ_muldet_init(m, ndets) < 0)
		return -1;

	struct iter_muldet_data id = { .i = 0, .muldet = m };
	if (data_multidet_foreach(fid, iter_muldet, &id) < 0)
		return -1;

	return 0;
}

void circ_prog_init(struct circ_prog *prog, size_t len)
{
	prog->i = 0;
	prog->len = len;
	prog->pc = 0;
}

void circ_prog_tick(struct circ_prog *prog)
{
	prog->i++;

	const unsigned pc = prog->i * 100 / prog->len;
	if (pc > prog->pc) {
		prog->pc = pc;
		log_info("Progress: %u%%", prog->pc);
	}
}

int circ_values_init(struct circ_values *vals, size_t len)
{
	_Complex double *z = malloc(sizeof(_Complex double) * len);
	if (!z)
		return -1;

	vals->z = z;
	vals->len = len;

	return 0;
}

void circ_values_free(struct circ_values *vals)
{
	free(vals->z);
}

int circ_init(struct circ *ct, const data_id fid, const size_t vals_len)
{
	memset(&ct->md, 0, sizeof ct->md);
	memset(&ct->cm, 0, sizeof ct->cm);

	if (circ_hamil_from_file(&ct->hm, fid) < 0) {
		log_error("circ_init: loading Hamiltonian failed");
		goto err_hamil_init;
	}
	log_debug("circ_init: Hamiltonian loaded (%u qubits, %zu terms)",
		ct->hm.qb, ct->hm.len);

	if (data_state_prep_kind(fid, &ct->stprep_kind) < 0) {
		log_error("circ_init: state_prep kind probe failed");
		goto err_stprep_kind;
	}

	switch (ct->stprep_kind) {
	case STPREP_MULTIDET:
		if (circ_muldet_from_file(&ct->md, fid) < 0) {
			log_error("circ_init: loading multidet state failed");
			goto err_stprep_load;
		}
		log_debug("circ_init: multidet state (%zu dets)", ct->md.len);
		break;
	case STPREP_COEFF_MATRIX:
		if (circ_coeff_init(&ct->cm, fid) < 0) {
			log_error("circ_init: coeff_matrix init failed");
			goto err_stprep_load;
		}
		log_debug("circ_init: coeff_matrix state (n_components=%zu)",
			ct->cm.n_components);
		break;
	}

	if (qreg_init(&ct->reg, ct->hm.qb) < 0) {
		log_error("circ_init: qreg_init failed");
		goto err_qreg_init;
	}
	if (circ_cache_init(ct->reg.qb_hi, ct->reg.qb_lo) < 0) {
		log_error("circ_init: cache_init failed");
		goto err_cache_init;
	}
	if (circ_values_init(&ct->vals, vals_len) < 0) {
		log_error("circ_init: values_init failed (len=%zu)", vals_len);
		goto err_vals_init;
	}

	return 0;

	circ_values_free(&ct->vals);
err_vals_init:
err_cache_init:
	qreg_free(&ct->reg);
err_qreg_init:
	switch (ct->stprep_kind) {
	case STPREP_MULTIDET:
		circ_muldet_free(&ct->md);
		break;
	case STPREP_COEFF_MATRIX:
		circ_coeff_free(&ct->cm);
		break;
	}
err_stprep_load:
err_stprep_kind:
	circ_hamil_free(&ct->hm);
err_hamil_init:
	return -1;
}

void circ_free(struct circ *ct)
{
	circ_values_free(&ct->vals);
	circ_hamil_free(&ct->hm);
	switch (ct->stprep_kind) {
	case STPREP_MULTIDET:
		circ_muldet_free(&ct->md);
		break;
	case STPREP_COEFF_MATRIX:
		circ_coeff_free(&ct->cm);
		break;
	}
	qreg_free(&ct->reg);
}

int circ_prepst(struct circ *ct)
{
	qreg_zero(&ct->reg);

	switch (ct->stprep_kind) {
	case STPREP_MULTIDET: {
		const struct circ_muldet *md = &ct->md;
		for (size_t i = 0; i < md->len; i++)
			qreg_setamp(&ct->reg, md->dets[i].idx, md->dets[i].cf);
		return 0;
	}
	case STPREP_COEFF_MATRIX:
		return state_prep_coeff_expand_all(&ct->reg, &ct->cm);
	}

	return -1;
}

static void circ_flush(struct paulis code_hi, const struct paulis *codes_lo,
	double *phis, size_t ncodes, void *data)
{
	struct qreg *reg = data;
	qreg_paulirot(reg, code_hi, codes_lo, phis, ncodes);
}

/*
 * circ_step_generic - apply one Trotter step of the
 * Hamiltonian evolution exp(i*omega*H).
 *
 * Each term is split into hi and lo Pauli codes.  If the
 * hi code matches the current cache group, the term is
 * appended (sharing the MPI exchange).  If the hi code
 * differs or the cache is full, the cache is flushed
 * (triggering MPI exchange + batched lo-rotations), then
 * the new term starts a fresh cache group.
 *
 * The Hamiltonian should be pre-sorted so that terms with
 * identical hi codes are contiguous — see circ_hamil_sort_lex.
 */
static int circ_step_generic(struct circ *ct, const struct circ_hamil *hm,
	const double omega, bool reverse)
{
	for (size_t i = 0; i < hm->len; i++) {
		size_t j = i;
		if (reverse)
			j = hm->len - i - 1;
		const double phi = omega * hm->terms[j].cf;
		const struct paulis code = hm->terms[j].op;

		if (circ_cache_insert(code, phi) == 0)
			continue;

		log_trace("paulirot, term: %zu, num_codes: %zu", i,
			circ_cache_len());
		circ_cache_flush(circ_flush, &ct->reg);
		if (circ_cache_insert(code, phi) < 0)
			return -1;
	}
	log_trace(
		"paulirot, last term group, num_codes: %zu", circ_cache_len());
	circ_cache_flush(circ_flush, &ct->reg);

	return 0;
}

inline int circ_step(
	struct circ *ct, const struct circ_hamil *hm, const double omega)
{
	return circ_step_generic(ct, hm, omega, false);
}

inline int circ_step_reverse(
	struct circ *ct, const struct circ_hamil *hm, const double omega)
{
	return circ_step_generic(ct, hm, omega, true);
}

/*
 * Inner product <trial | evolved> for the coeff_matrix path.
 *
 * Walks the same Slater-Condon outer product used at expand
 * time, summing conj(trial(idx)) * evolved(idx) over the
 * amplitudes owned by this rank, then MPI-reduces.  The
 * walk costs O(M_alpha * M_beta) per measurement, on the
 * same order as expansion itself.
 */
static _Complex double measure_coeff_block(struct qreg *reg,
	const uint32_t n_sites, const uint32_t n_alpha, const uint32_t n_beta,
	const double *C_alpha, const double *C_beta, const double weight,
	const int tapered)
{
	return state_prep_coeff_inner(reg, n_sites, n_alpha, n_beta, C_alpha,
		C_beta, weight, tapered);
}

static _Complex double measure_coeff(struct circ *ct)
{
	const struct circ_coeff *cm = &ct->cm;
	_Complex double pr = 0.0;

	if (cm->n_components == 0) {
		pr = measure_coeff_block(&ct->reg, cm->n_sites, cm->n_alpha,
			cm->n_beta, cm->C_alpha,
			cm->closed_shell ? NULL : cm->C_beta, 1.0,
			cm->tapered);
	} else {
		for (size_t k = 0; k < cm->n_components; k++) {
			const struct circ_coeff_block *b = &cm->blocks[k];
			pr += measure_coeff_block(&ct->reg, cm->n_sites,
				cm->n_alpha, cm->n_beta, b->C_alpha,
				cm->closed_shell ? NULL : b->C_beta, b->cf,
				cm->tapered);
		}
	}

	return pr;
}

_Complex double circ_measure(struct circ *ct)
{
	switch (ct->stprep_kind) {
	case STPREP_MULTIDET: {
		const struct circ_muldet *md = &ct->md;
		_Complex double pr = 0.0;
		for (size_t i = 0; i < md->len; i++) {
			_Complex double a;
			qreg_getamp(&ct->reg, md->dets[i].idx, &a);
			pr += a * conj(md->dets[i].cf);
		}
		return pr;
	}
	case STPREP_COEFF_MATRIX:
		return measure_coeff(ct);
	}

	return 0.0;
}
