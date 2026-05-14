/*
 * Slater-Condon expansion of a coefficient-matrix trial
 * state.  See include/phase2/state_prep_coeff.h for the
 * public API contract; this file holds the implementation.
 */

#define LOG_SUBSYS "prep"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "combinations.h"
#include "det_small.h"
#include "log.h"
#include "phase2/data.h"
#include "phase2/qreg.h"
#include "phase2/state_prep_coeff.h"

#include "qreg.h"

#define SPARSITY_PRUNE (1.0e-12)

static inline uint64_t occ_pair_to_idx(const uint32_t *occ_a,
	const uint32_t na, const uint32_t *occ_b, const uint32_t nb,
	const uint32_t n_sites)
{
	uint64_t idx = 0;
	for (uint32_t i = 0; i < na; i++)
		idx |= (UINT64_C(1) << occ_a[i]);
	for (uint32_t i = 0; i < nb; i++)
		idx |= (UINT64_C(1) << (n_sites + occ_b[i]));
	return idx;
}

static inline uint64_t drop_two_bits(const uint64_t idx, const uint32_t n_sites)
{
	const uint64_t lo = (idx >> 1) & ((UINT64_C(1) << (n_sites - 1)) - 1);
	const uint64_t hi = idx >> (n_sites + 1);
	return lo | (hi << (n_sites - 1));
}

static uint64_t binomial(const uint32_t n, const uint32_t k)
{
	if (k > n)
		return 0;
	uint32_t r = k;
	if (r > n - r)
		r = n - r;
	uint64_t num = 1, den = 1;
	for (uint32_t i = 1; i <= r; i++) {
		num *= (uint64_t)(n - r + i);
		den *= (uint64_t)i;
	}
	return num / den;
}

/*
 * Enumerate all k-subsets in lex order, fill out[]
 * (size M*k) with the indices, and dets[] (size M) with
 * det(C[occ, :]) for each subset.
 *
 * C is row-major (n_sites, k).  Tuple `occ` selects rows.
 */
static int precompute_dets(uint32_t n_sites, uint32_t k, const double *C,
	uint32_t *out_tuples, double *out_dets, size_t M)
{
	struct combo it;
	combinations_init(&it, n_sites, k);

	uint32_t buf[COMBINATIONS_MAX_K];
	double sub[DET_SMALL_MAX_N * DET_SMALL_MAX_N];

	size_t i = 0;
	while (combinations_next(&it, buf) == 0) {
		if (i >= M)
			return -1;
		if (k == 0) {
			out_dets[i] = 1.0;
		} else {
			for (uint32_t r = 0; r < k; r++) {
				const uint32_t row = buf[r];
				for (uint32_t c = 0; c < k; c++)
					sub[r * k + c] = C[row * k + c];
			}
			out_dets[i] = det_small(sub, k);
		}
		for (uint32_t r = 0; r < k; r++)
			out_tuples[i * (k ? k : 1) + r] = buf[r];
		i++;
	}

	if (i != M)
		return -1;
	return 0;
}

int state_prep_coeff_expand(struct qreg *reg, const uint32_t n_sites,
	const uint32_t n_alpha, const uint32_t n_beta, const double *C_alpha,
	const double *C_beta, const double weight, const int tapered,
	const int accumulate)
{
	if (n_alpha > n_sites || n_beta > n_sites)
		return -1;
	if (n_alpha > DET_SMALL_MAX_N || n_beta > DET_SMALL_MAX_N)
		return -1;
	if (n_alpha > COMBINATIONS_MAX_K || n_beta > COMBINATIONS_MAX_K)
		return -1;

	const double *Cb = C_beta ? C_beta : C_alpha;

	const size_t Ma = (size_t)binomial(n_sites, n_alpha);
	const size_t Mb = (size_t)binomial(n_sites, n_beta);

	const uint32_t ka = n_alpha ? n_alpha : 1;
	const uint32_t kb = n_beta ? n_beta : 1;

	int rt = -1;

	uint32_t *tup_a = malloc(sizeof(uint32_t) * Ma * ka);
	if (!tup_a)
		return -1;
	uint32_t *tup_b = malloc(sizeof(uint32_t) * Mb * kb);
	if (!tup_b)
		goto err_tb;
	double *det_a = malloc(sizeof(double) * Ma);
	if (!det_a)
		goto err_da;
	double *det_b = malloc(sizeof(double) * Mb);
	if (!det_b)
		goto err_db;

	if (precompute_dets(n_sites, n_alpha, C_alpha, tup_a, det_a, Ma) < 0)
		goto err_pre;
	if (precompute_dets(n_sites, n_beta, Cb, tup_b, det_b, Mb) < 0)
		goto err_pre;

	const int my_rank = reg->wd.rank;

	for (size_t i = 0; i < Ma; i++) {
		const double da = det_a[i];
		const uint32_t *oa = &tup_a[i * ka];
		for (size_t j = 0; j < Mb; j++) {
			const double cf = weight * da * det_b[j];
			if (fabs(cf) < SPARSITY_PRUNE)
				continue;
			const uint32_t *ob = &tup_b[j * kb];
			const uint64_t full = occ_pair_to_idx(
				oa, n_alpha, ob, n_beta, n_sites);
			const uint64_t idx =
				tapered ? drop_two_bits(full, n_sites) : full;

			const uint64_t owner = qreg_getihi(reg, idx);
			if (owner != (uint64_t)my_rank)
				continue;
			const uint64_t i_lo = qreg_getilo(reg, idx);
			if (accumulate)
				reg->amp[i_lo] += cf;
			else
				reg->amp[i_lo] = cf;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	rt = 0;
err_pre:
	free(det_b);
err_db:
	free(det_a);
err_da:
	free(tup_b);
err_tb:
	free(tup_a);
	return rt;
}

_Complex double state_prep_coeff_inner(struct qreg *reg,
	const uint32_t n_sites, const uint32_t n_alpha, const uint32_t n_beta,
	const double *C_alpha, const double *C_beta, const double weight,
	const int tapered)
{
	if (n_alpha > n_sites || n_beta > n_sites)
		return 0.0;
	if (n_alpha > DET_SMALL_MAX_N || n_beta > DET_SMALL_MAX_N)
		return 0.0;
	if (n_alpha > COMBINATIONS_MAX_K || n_beta > COMBINATIONS_MAX_K)
		return 0.0;

	const double *Cb = C_beta ? C_beta : C_alpha;
	const size_t Ma = (size_t)binomial(n_sites, n_alpha);
	const size_t Mb = (size_t)binomial(n_sites, n_beta);
	const uint32_t ka = n_alpha ? n_alpha : 1;
	const uint32_t kb = n_beta ? n_beta : 1;

	uint32_t *tup_a = malloc(sizeof(uint32_t) * Ma * ka);
	uint32_t *tup_b = malloc(sizeof(uint32_t) * Mb * kb);
	double *det_a = malloc(sizeof(double) * Ma);
	double *det_b = malloc(sizeof(double) * Mb);
	if (!tup_a || !tup_b || !det_a || !det_b) {
		free(tup_a);
		free(tup_b);
		free(det_a);
		free(det_b);
		return 0.0;
	}
	if (precompute_dets(n_sites, n_alpha, C_alpha, tup_a, det_a, Ma) < 0
		|| precompute_dets(n_sites, n_beta, Cb, tup_b, det_b, Mb)
			< 0) {
		free(tup_a);
		free(tup_b);
		free(det_a);
		free(det_b);
		return 0.0;
	}

	const int my_rank = reg->wd.rank;

	double acc_r = 0.0, acc_i = 0.0;
	for (size_t i = 0; i < Ma; i++) {
		const double da = det_a[i];
		const uint32_t *oa = &tup_a[i * ka];
		for (size_t j = 0; j < Mb; j++) {
			const double cf = weight * da * det_b[j];
			if (fabs(cf) < SPARSITY_PRUNE)
				continue;
			const uint32_t *ob = &tup_b[j * kb];
			const uint64_t full = occ_pair_to_idx(
				oa, n_alpha, ob, n_beta, n_sites);
			const uint64_t idx =
				tapered ? drop_two_bits(full, n_sites) : full;
			const uint64_t owner = qreg_getihi(reg, idx);
			if (owner != (uint64_t)my_rank)
				continue;
			const uint64_t i_lo = qreg_getilo(reg, idx);
			const _Complex double a = reg->amp[i_lo];
			acc_r += cf * creal(a);
			acc_i += cf * cimag(a);
		}
	}

	double partial[2] = { acc_r, acc_i };
	double total[2] = { 0.0, 0.0 };
	MPI_Allreduce(partial, total, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	free(tup_a);
	free(tup_b);
	free(det_a);
	free(det_b);

	return CMPLX(total[0], total[1]);
}

int state_prep_coeff_expand_all(struct qreg *reg, const struct circ_coeff *cm)
{
	log_debug("expand_all: n_sites=%u n_alpha=%u n_beta=%u"
		  " closed_shell=%d tapered=%d n_components=%zu",
		cm->n_sites, cm->n_alpha, cm->n_beta, cm->closed_shell,
		cm->tapered, cm->n_components);

	if (cm->n_components == 0) {
		if (state_prep_coeff_expand(reg, cm->n_sites, cm->n_alpha,
			    cm->n_beta, cm->C_alpha,
			    cm->closed_shell ? NULL : cm->C_beta, 1.0,
			    cm->tapered, 0) < 0) {
			log_error("expand_all: single-block expand failed");
			return -1;
		}
		return 0;
	}

	qreg_zero(reg);
	for (size_t k = 0; k < cm->n_components; k++) {
		const struct circ_coeff_block *b = &cm->blocks[k];
		log_trace("expand_all: block %zu/%zu cf=%.6f", k + 1,
			cm->n_components, b->cf);
		if (state_prep_coeff_expand(reg, cm->n_sites, cm->n_alpha,
			    cm->n_beta, b->C_alpha,
			    cm->closed_shell ? NULL : b->C_beta, b->cf,
			    cm->tapered, 1) < 0) {
			log_error("expand_all: block %zu expand failed", k);
			return -1;
		}
	}
	return 0;
}

int circ_coeff_init(struct circ_coeff *cm, const data_id fid)
{
	memset(cm, 0, sizeof *cm);

	if (data_coeff_matrix_getnums(fid, &cm->n_qubits, &cm->n_sites,
		    &cm->n_alpha, &cm->n_beta, &cm->closed_shell,
		    &cm->tapered) < 0)
		return -1;

	size_t n_comp = 0;
	if (data_coeff_matrix_csf_count(fid, &n_comp) < 0)
		return -1;

	const size_t sz_a = (size_t)cm->n_sites * cm->n_alpha;
	const size_t sz_b = (size_t)cm->n_sites * cm->n_beta;

	cm->C_alpha = malloc(sizeof(double) * (sz_a ? sz_a : 1));
	if (!cm->C_alpha)
		goto err_ca;
	if (!cm->closed_shell) {
		cm->C_beta = malloc(sizeof(double) * (sz_b ? sz_b : 1));
		if (!cm->C_beta)
			goto err_cb;
	}

	if (n_comp == 0) {
		if (data_coeff_matrix_read(fid, cm->C_alpha, cm->C_beta) < 0)
			goto err_read;
		cm->n_components = 0;
		return 0;
	}

	cm->n_components = n_comp;
	cm->blocks = calloc(n_comp, sizeof *cm->blocks);
	if (!cm->blocks)
		goto err_blocks;

	for (size_t k = 0; k < n_comp; k++) {
		cm->blocks[k].C_alpha =
			malloc(sizeof(double) * (sz_a ? sz_a : 1));
		if (!cm->blocks[k].C_alpha)
			goto err_inner;
		if (!cm->closed_shell) {
			cm->blocks[k].C_beta =
				malloc(sizeof(double) * (sz_b ? sz_b : 1));
			if (!cm->blocks[k].C_beta)
				goto err_inner;
		}
		if (data_coeff_matrix_csf_read(fid, k, &cm->blocks[k].cf,
			    cm->blocks[k].C_alpha,
			    cm->closed_shell ? NULL
					     : cm->blocks[k].C_beta) < 0)
			goto err_inner;
	}
	return 0;

err_inner:
	for (size_t k = 0; k < n_comp; k++) {
		free(cm->blocks[k].C_alpha);
		free(cm->blocks[k].C_beta);
	}
	free(cm->blocks);
	cm->blocks = NULL;
err_blocks:
err_read:
	if (cm->C_beta) {
		free(cm->C_beta);
		cm->C_beta = NULL;
	}
err_cb:
	free(cm->C_alpha);
	cm->C_alpha = NULL;
err_ca:
	return -1;
}

void circ_coeff_free(struct circ_coeff *cm)
{
	if (!cm)
		return;
	if (cm->blocks) {
		for (size_t k = 0; k < cm->n_components; k++) {
			free(cm->blocks[k].C_alpha);
			free(cm->blocks[k].C_beta);
		}
		free(cm->blocks);
		cm->blocks = NULL;
	}
	free(cm->C_alpha);
	cm->C_alpha = NULL;
	free(cm->C_beta);
	cm->C_beta = NULL;
}
