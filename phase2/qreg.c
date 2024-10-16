#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"

#define MAX_COUNT (1 << 29)

static int qreg_backend_init(struct qreg *reg);
static int qreg_backend_destroy(struct qreg *reg);

uint64_t qreg_getilo(const struct qreg *reg, uint64_t i)
{
	const uint64_t mask_lo = (UINT64_C(1) << reg->nqb_lo) - 1;

	return i & mask_lo;
}

uint64_t qreg_getihi(const struct qreg *reg, uint64_t i)
{
	const uint64_t mask_hi = (UINT64_C(1) << reg->nqb_hi) - 1;

	return (i >> reg->nqb_lo) & mask_hi;
}

int qreg_init(struct qreg *reg, const uint32_t nqb)
{
	if (world_info(&reg->wd) != WORLD_READY)
		return -1;

	uint32_t nqb_hi = 0, nrk = reg->wd.size;
	while (nrk >>= 1) /* nqb_hi = log2(nrk) */
		nqb_hi++;
	if (nqb_hi >= nqb)
		return -1;
	const uint32_t nqb_lo = nqb - nqb_hi;
	const uint64_t namp = UINT64_C(1) << nqb_lo;

	const int msg_count = namp < MAX_COUNT ? namp : MAX_COUNT;
	const size_t nreqs = namp / msg_count;

	MPI_Request *const reqs = malloc(sizeof *reqs * nreqs * 2);
	if (reqs == nullptr)
		goto err_reqs_alloc;
	_Complex double *const amp = malloc(sizeof *amp * namp * 2);
	if (amp == nullptr)
		goto err_amp_alloc;

	reg->nqb_lo = nqb_lo;
	reg->nqb_hi = nqb_hi;
	reg->amp = amp;
	reg->buf = amp + namp;
	reg->namp = namp;
	reg->msg_count = msg_count;
	reg->reqs_snd = reqs;
	reg->reqs_rcv = reqs + nreqs;
	reg->nreqs = nreqs;

	if (qreg_backend_init(reg) < 0)
		goto err_backend_init;

	return 0;

err_backend_init:
	free(amp);
err_amp_alloc:
	free(reqs);
err_reqs_alloc:
	return -1;
}

void qreg_destroy(struct qreg *reg)
{
	qreg_backend_destroy(reg);

	if (reg->amp != nullptr)
		free(reg->amp);
	if (reg->reqs_snd != nullptr)
		free(reg->reqs_snd);
}

#if PHASE2_BACKEND == 0 /* qreg */

static __inline__ int qreg_backend_init(struct qreg *reg)
{
	(void)reg;

	return 0;
}

static __inline__ int qreg_backend_destroy(struct qreg *reg)
{
	(void)reg;

	return 0;
}

#elif PHASE2_BACKEND == 1 /* QuEST */
int qreg_quest_init(struct qreg *reg);
int qreg_quest_destroy(struct qreg *reg);

static __inline__ int qreg_backend_init(struct qreg *reg)
{
	return qreg_quest_init(reg);
}

static __inline__ int qreg_backend_destroy(struct qreg *reg)
{
	return qreg_quest_destroy(reg);
}

#elif PHASE2_BACKEND == 2 /* CUDA */
int qreg_cuda_init(struct qreg *reg);
int qreg_cuda_destroy(struct qreg *reg);

static __inline__ int qreg_backend_init(struct qreg *reg)
{
	return qreg_cuda_init(reg);
}

static __inline__ int qreg_backend_destroy(struct qreg *reg)
{
	return qreg_cuda_destroy(reg);
}
#endif /* PHASE2_BACKEND */
