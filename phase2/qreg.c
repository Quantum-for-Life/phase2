/*
 * Backend-neutral qreg lifecycle: amplitude buffer
 * allocation, hi/lo index split (qb_hi = log2(MPI
 * ranks), qb_lo = nqb - qb_hi), MPI request scratch.
 * The two backends (qreg_qreg.c CPU, qreg_cuda.c GPU)
 * supply the operator kernels via the private hooks
 * declared in phase2/qreg.h.
 */

#define LOG_SUBSYS "qreg"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "log.h"
#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"

#include "qreg.h"

/* Per-MPI-message element count cap.  Each cdouble
 * is sent as 2 * MPI_DOUBLE; MPI_Isend's `count`
 * is int, so 2 * MAX_COUNT must fit in int.  At
 * 2^29 cdoubles (= 8 GiB / send) the per-chunk MPI
 * overhead is negligible and the doubles count
 * 2^30 fits the int boundary with headroom. */
#define MAX_COUNT (1 << 29)

uint64_t qreg_getilo(const struct qreg *reg, uint64_t i)
{
	const uint64_t mask_lo = (UINT64_C(1) << reg->qb_lo) - 1;

	return i & mask_lo;
}

uint64_t qreg_getihi(const struct qreg *reg, uint64_t i)
{
	const uint64_t mask_hi = (UINT64_C(1) << reg->qb_hi) - 1;

	return (i >> reg->qb_lo) & mask_hi;
}

int qreg_init(struct qreg *reg, const uint32_t nqb)
{
	if (world_info(&reg->wd) != WORLD_READY) {
		log_error("qreg_init: world not ready");
		return -1;
	}

	uint32_t nqb_hi = 0, nrk = reg->wd.size;
	while (nrk >>= 1) /* nqb_hi = log2(nrk) */
		nqb_hi++;
	if (nqb_hi >= nqb) {
		log_error("qreg_init: nqb=%u must exceed log2(ranks)=%u",
			nqb, nqb_hi);
		return -1;
	}
	const uint32_t nqb_lo = nqb - nqb_hi;
	const uint64_t namp = UINT64_C(1) << nqb_lo;

	const int msg_count = namp < MAX_COUNT ? namp : MAX_COUNT;
	const size_t nreqs = namp / msg_count;

	MPI_Request *const reqs = malloc(sizeof *reqs * nreqs * 2);
	if (!reqs) {
		log_error("qreg_init: alloc reqs (%zu bytes)",
			sizeof *reqs * nreqs * 2);
		goto err_reqs_alloc;
	}
	_Complex double *const amp = malloc(sizeof *amp * namp * 2);
	if (!amp) {
		log_error("qreg_init: alloc amp (%zu bytes)",
			sizeof *amp * namp * 2);
		goto err_amp_alloc;
	}

	reg->qb_lo = nqb_lo;
	reg->qb_hi = nqb_hi;
	reg->amp = amp;
	reg->buf = amp + namp;
	reg->namp = namp;
	reg->msg_count = msg_count;
	reg->reqs_snd = reqs;
	reg->reqs_rcv = reqs + nreqs;
	reg->nreqs = nreqs;

	if (qreg_backend_init(reg) < 0) {
		log_error("qreg_init: backend init failed");
		goto err_backend_init;
	}

	log_debug("qreg_init: nqb=%u qb_lo=%u qb_hi=%u namp=%llu"
		  " msg_count=%d nreqs=%zu",
		nqb, nqb_lo, nqb_hi, (unsigned long long)namp, msg_count,
		nreqs);
	return 0;

err_backend_init:
	free(amp);
err_amp_alloc:
	free(reqs);
err_reqs_alloc:
	return -1;
}

void qreg_free(struct qreg *reg)
{
	qreg_backend_free(reg);

	if (reg->amp != nullptr)
		free(reg->amp);
	if (reg->reqs_snd != nullptr)
		free(reg->reqs_snd);
}

/* Resolve partner rank from code_hi, post the exchange,
 * compute the phase bm carried from the partner, wait. */
static void qreg_paulirot_hi(struct qreg *reg, struct paulis code_hi,
	_Complex double *bm)
{
	paulis_shr(&code_hi, reg->qb_lo);
	const uint64_t rnk_rem = paulis_effect(code_hi, reg->wd.rank, nullptr);

	qreg_backend_exch_init(reg, rnk_rem);
	paulis_effect(code_hi, rnk_rem, bm);
	qreg_backend_exch_waitall(reg);
}

void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *phis,
	const size_t ncodes)
{
	_Complex double bm = 1.0;
	qreg_paulirot_hi(reg, code_hi, &bm);
	qreg_backend_paulirot_lo(reg, codes_lo, phis, ncodes, bm);
}
