#include "c23_compat.h"
#include <math.h>
#include <stdlib.h>

#include "mpi.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"

#include "qreg.h"

#define MAX_COUNT (1UL << 29)

uint64_t qreg_getilo(const struct qreg *reg, const uint64_t i)
{
	const uint64_t mask_lo = (UINT64_C(1) << reg->qb_lo) - 1;

	return i & mask_lo;
}

uint64_t qreg_getihi(const struct qreg *reg, const uint64_t i)
{
	const uint64_t mask_hi = (UINT64_C(1) << reg->qb_hi) - 1;

	return (i >> reg->qb_lo) & mask_hi;
}

int qreg_init(struct qreg *reg, const uint32_t qb)
{
	if (world_info(&reg->wd) != WORLD_READY)
		return -1;

	uint32_t qb_hi = 0, nrk = reg->wd.size;
	while (nrk >>= 1) /* qb_hi = log2(nrk) */
		qb_hi++;
	if (qb_hi >= qb)
		return -1;
	const uint32_t qb_lo = qb - qb_hi;
	const uint64_t namp = UINT64_C(1) << qb_lo;

	const int msg_count = namp < MAX_COUNT ? namp : MAX_COUNT;
	const size_t nreqs = namp / msg_count;

	MPI_Request *const reqs = malloc(sizeof *reqs * nreqs * 2);
	if (!reqs)
		goto err_reqs_alloc;
	_Complex double *const amp = malloc(sizeof *amp * namp * 2);
	if (!amp)
		goto err_amp_alloc;

	reg->qb_lo = qb_lo;
	reg->qb_hi = qb_hi;
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
