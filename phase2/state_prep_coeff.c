/*
 * Slater-Condon expansion of a coefficient-matrix trial
 * state, plus the matching inner-product readout.  See
 * include/phase2/state_prep_coeff.h for the public API
 * and the scratch contract.
 */

#define LOG_SUBSYS "prep"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mpi.h"

#include "combinations.h"
#include "det_small.h"
#include "log.h"
#include "phase2/circ.h"
#include "phase2/qreg.h"
#include "phase2/state_prep_coeff.h"

#include "qreg.h"

/* Coefficient magnitudes below this are skipped at
 * expand / inner time (sparsity prune). */
#define SPARSITY_PRUNE (1.0e-12)

/* Pack alpha occupations into the low n_sites bits and
 * beta occupations into the next n_sites bits. */
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

/* Tapered-state index transform: drop bit 0 (the alpha
 * symmetry generator) and bit n_sites (the beta one),
 * compacting the remaining bits into 2*(n_sites-1). */
static inline uint64_t drop_two_bits(const uint64_t idx, const uint32_t n_sites)
{
	const uint64_t lo = (idx >> 1) & ((UINT64_C(1) << (n_sites - 1)) - 1);
	const uint64_t hi = idx >> (n_sites + 1);
	return lo | (hi << (n_sites - 1));
}

/* C(n, k); 0 when k > n. */
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

/* Fill out[] (size M*k_eff) with the k-subset tuples of
 * (n_sites choose k) in lex order.  k_eff is max(k, 1)
 * so the M=1 / k=0 degenerate case still occupies one
 * row.  Returns 0 on success, -1 if the iterator yields
 * the wrong number of tuples. */
static int precompute_tuples(uint32_t n_sites, uint32_t k,
	uint32_t *out_tuples, size_t M)
{
	struct combo it;
	combinations_init(&it, n_sites, k);

	uint32_t buf[COMBINATIONS_MAX_K];
	const uint32_t k_eff = k ? k : 1;
	size_t i = 0;
	while (combinations_next(&it, buf) == 0) {
		if (i >= M)
			return -1;
		for (uint32_t r = 0; r < k; r++)
			out_tuples[i * k_eff + r] = buf[r];
		i++;
	}
	return i == M ? 0 : -1;
}

/* Walk pre-filled tuples and compute det(C[occ, :]) for
 * each.  C is row-major (n_sites, k).  k=0 degenerates
 * to det = 1. */
static void compute_dets(uint32_t k, const uint32_t *tuples,
	const double *C, double *out_dets, size_t M)
{
	if (k == 0) {
		for (size_t i = 0; i < M; i++)
			out_dets[i] = 1.0;
		return;
	}
	double sub[DET_SMALL_MAX_N * DET_SMALL_MAX_N];
	for (size_t i = 0; i < M; i++) {
		const uint32_t *occ = &tuples[i * k];
		for (uint32_t r = 0; r < k; r++) {
			const uint32_t row = occ[r];
			for (uint32_t c = 0; c < k; c++)
				sub[r * k + c] = C[row * k + c];
		}
		out_dets[i] = det_small(sub, k);
	}
}


int state_prep_coeff_scratch_init(struct state_prep_coeff_scratch *sc,
	uint32_t n_sites, uint32_t n_alpha, uint32_t n_beta)
{
	if (n_alpha > n_sites || n_beta > n_sites) {
		log_error("scratch_init: n_alpha=%u or n_beta=%u > n_sites=%u",
			n_alpha, n_beta, n_sites);
		return -1;
	}
	if (n_alpha > DET_SMALL_MAX_N || n_beta > DET_SMALL_MAX_N) {
		log_error("scratch_init: n > DET_SMALL_MAX_N=%d",
			DET_SMALL_MAX_N);
		return -1;
	}
	if (n_alpha > COMBINATIONS_MAX_K || n_beta > COMBINATIONS_MAX_K) {
		log_error("scratch_init: n > COMBINATIONS_MAX_K=%d",
			COMBINATIONS_MAX_K);
		return -1;
	}

	memset(sc, 0, sizeof *sc);
	sc->n_sites = n_sites;
	sc->n_alpha = n_alpha;
	sc->n_beta = n_beta;
	sc->Ma = (size_t)binomial(n_sites, n_alpha);
	sc->Mb = (size_t)binomial(n_sites, n_beta);

	const uint32_t ka = n_alpha ? n_alpha : 1;
	const uint32_t kb = n_beta ? n_beta : 1;

	sc->tup_a = malloc(sizeof *sc->tup_a * sc->Ma * ka);
	sc->tup_b = malloc(sizeof *sc->tup_b * sc->Mb * kb);
	sc->det_a = malloc(sizeof *sc->det_a * sc->Ma);
	sc->det_b = malloc(sizeof *sc->det_b * sc->Mb);
	if (!sc->tup_a || !sc->tup_b || !sc->det_a || !sc->det_b) {
		log_error("scratch_init: alloc failed");
		goto err_alloc;
	}

	if (precompute_tuples(n_sites, n_alpha, sc->tup_a, sc->Ma) < 0)
		goto err_alloc;
	if (precompute_tuples(n_sites, n_beta, sc->tup_b, sc->Mb) < 0)
		goto err_alloc;

	return 0;

err_alloc:
	state_prep_coeff_scratch_free(sc);
	return -1;
}

void state_prep_coeff_scratch_free(struct state_prep_coeff_scratch *sc)
{
	free(sc->tup_a);
	free(sc->tup_b);
	free(sc->det_a);
	free(sc->det_b);
	memset(sc, 0, sizeof *sc);
}


int state_prep_coeff_expand(struct qreg *reg,
	struct state_prep_coeff_scratch *sc,
	const double *C_alpha, const double *C_beta,
	const double weight, const int tapered, const int accumulate)
{
	const double *Cb = C_beta ? C_beta : C_alpha;

	compute_dets(sc->n_alpha, sc->tup_a, C_alpha, sc->det_a, sc->Ma);
	compute_dets(sc->n_beta, sc->tup_b, Cb, sc->det_b, sc->Mb);

	const int my_rank = reg->wd.rank;
	const uint32_t ka = sc->n_alpha ? sc->n_alpha : 1;
	const uint32_t kb = sc->n_beta ? sc->n_beta : 1;

	for (size_t i = 0; i < sc->Ma; i++) {
		const double da = sc->det_a[i];
		const uint32_t *oa = &sc->tup_a[i * ka];
		for (size_t j = 0; j < sc->Mb; j++) {
			const double cf = weight * da * sc->det_b[j];
			if (fabs(cf) < SPARSITY_PRUNE)
				continue;
			const uint32_t *ob = &sc->tup_b[j * kb];
			const uint64_t full = occ_pair_to_idx(
				oa, sc->n_alpha, ob, sc->n_beta, sc->n_sites);
			const uint64_t idx = tapered
				? drop_two_bits(full, sc->n_sites)
				: full;

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
	return 0;
}

int state_prep_coeff_inner(struct qreg *reg,
	struct state_prep_coeff_scratch *sc,
	const double *C_alpha, const double *C_beta,
	const double weight, const int tapered, _Complex double *out)
{
	const double *Cb = C_beta ? C_beta : C_alpha;

	/* On the CUDA backend the canonical state lives in
	 * cu->damp; reg->amp may be stale after any kernel
	 * run (circuit evolution, paulirot, etc.).  Pull the
	 * device slab back to the host before reading.  CPU
	 * backend: no-op + barrier. */
	qreg_sync_device_to_host(reg);

	compute_dets(sc->n_alpha, sc->tup_a, C_alpha, sc->det_a, sc->Ma);
	compute_dets(sc->n_beta, sc->tup_b, Cb, sc->det_b, sc->Mb);

	const int my_rank = reg->wd.rank;
	const uint32_t ka = sc->n_alpha ? sc->n_alpha : 1;
	const uint32_t kb = sc->n_beta ? sc->n_beta : 1;

	/* The outer walk is a pure map-reduce: every body
	 * local is declared in-loop (hence thread-private),
	 * the reads (det_a / det_b, tup_a / tup_b, reg->amp)
	 * are const for the duration, and the only carried
	 * state is the acc_r / acc_i sum -> an OpenMP
	 * reduction.  Builds without -fopenmp ignore the
	 * pragma and run the serial loop unchanged.  The
	 * reduction reorders the cross-thread partial sums,
	 * so the result can differ from the serial order in
	 * the last ULPs; this is an O(Ma * Mb) measurement
	 * accumulator whose terms are already summed in an
	 * MPI_Allreduce below (itself order-unspecified
	 * across ranks), so floating-point associativity is
	 * not relied upon here.  At OMP_NUM_THREADS=1 the
	 * order is identical to the serial loop. */
	double acc_r = 0.0, acc_i = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : acc_r, acc_i)
#endif
	for (size_t i = 0; i < sc->Ma; i++) {
		const double da = sc->det_a[i];
		const uint32_t *oa = &sc->tup_a[i * ka];
		for (size_t j = 0; j < sc->Mb; j++) {
			const double cf = weight * da * sc->det_b[j];
			if (fabs(cf) < SPARSITY_PRUNE)
				continue;
			const uint32_t *ob = &sc->tup_b[j * kb];
			const uint64_t full = occ_pair_to_idx(
				oa, sc->n_alpha, ob, sc->n_beta, sc->n_sites);
			const uint64_t idx = tapered
				? drop_two_bits(full, sc->n_sites)
				: full;
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
	*out = CMPLX(total[0], total[1]);
	return 0;
}

int state_prep_coeff_expand_all(struct qreg *reg,
	struct state_prep_coeff_scratch *sc,
	const struct data_coeff_matrix *cm)
{
	log_debug("expand_all: n_sites=%u n_alpha=%u n_beta=%u"
		  " closed_shell=%d tapered=%d n_components=%zu",
		cm->n_sites, cm->n_alpha, cm->n_beta, cm->closed_shell,
		cm->tapered, cm->n_components);

	qreg_zero(reg);

	/* qreg_zero clears the canonical state (reg->amp on
	 * CPU, cu->damp on CUDA).  state_prep_coeff_expand
	 * accumulates into reg->amp directly for
	 * performance, so the host buffer must also start
	 * from a known-zero state on the CUDA backend (where
	 * qreg_zero touches the device buffer only). */
	memset(reg->amp, 0,
		reg->namp * sizeof(_Complex double));

	if (cm->n_components == 0) {
		if (state_prep_coeff_expand(reg, sc, cm->C_alpha,
			    cm->closed_shell ? NULL : cm->C_beta, 1.0,
			    cm->tapered, 0) < 0) {
			log_error("expand_all: single-block expand failed");
			return -1;
		}
		qreg_sync_host_to_device(reg);
		return 0;
	}

	for (size_t k = 0; k < cm->n_components; k++) {
		const struct data_coeff_block *b = &cm->blocks[k];
		log_trace("expand_all: block %zu/%zu cf=%.6f", k + 1,
			cm->n_components, b->cf);
		if (state_prep_coeff_expand(reg, sc, b->C_alpha,
			    cm->closed_shell ? NULL : b->C_beta, b->cf,
			    cm->tapered, 1) < 0) {
			log_error("expand_all: block %zu expand failed", k);
			return -1;
		}
	}

	/* Push the fully-accumulated host buffer to the
	 * canonical state.  One bulk transfer per pak;
	 * preserves the tight host-side accumulation loop's
	 * performance. */
	qreg_sync_host_to_device(reg);
	return 0;
}
