#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"

#include "qreg.h"

typedef _Complex double c64;

int qreg_backend_init(struct qreg *reg)
{
	reg->data = nullptr;

	return 0;
}

void qreg_backend_destroy(struct qreg *reg)
{
	(void)reg;
}

void qreg_getamp(struct qreg *reg, const uint64_t i, c64 *z)
{
	const uint64_t i_lo = qreg_getilo(reg, i);
	const uint64_t rank = qreg_getihi(reg, i);

	if (reg->wd.rank == (int)rank)
		*z = reg->amp[i_lo];
	MPI_Bcast(z, 2, MPI_DOUBLE, rank, MPI_COMM_WORLD);
}

void qreg_setamp(struct qreg *reg, const uint64_t i, c64 z)
{
	const uint64_t i_lo = qreg_getilo(reg, i);
	const uint64_t rank = qreg_getihi(reg, i);

	if (reg->wd.rank == (int)rank)
		reg->amp[i_lo] = z;
	MPI_Barrier(MPI_COMM_WORLD);
}

void qreg_zero(struct qreg *reg)
{
	memset(reg->amp, 0, reg->namp * sizeof(c64));
}

static void exch_init(struct qreg *reg, const int rnk_rem)
{
	const int nr = reg->nreqs;

	for (int i = 0; i < nr; i++) {
		const size_t offset = i * reg->msg_count;

		MPI_Isend(reg->amp + offset, reg->msg_count * 2, MPI_DOUBLE,
			rnk_rem, i, MPI_COMM_WORLD, reg->reqs_snd + i);
		MPI_Irecv(reg->buf + offset, reg->msg_count * 2, MPI_DOUBLE,
			rnk_rem, i, MPI_COMM_WORLD, reg->reqs_rcv + i);
	}
}

static void exch_waitall(struct qreg *reg)
{
	MPI_Waitall(reg->nreqs, reg->reqs_snd, MPI_STATUSES_IGNORE);
	MPI_Waitall(reg->nreqs, reg->reqs_rcv, MPI_STATUSES_IGNORE);
}

/* These kernels can be easily ported to CUDA. */
static __inline__ void kernel_mix(
	const size_t i, c64 *restrict a, c64 *restrict b, const c64 bm)
{
	b[i] *= bm;

	const c64 x = a[i];
	const c64 y = b[i];
	a[i] = (x + y) / 2.0;
	b[i] = (x - y) / 2.0;
}

static __inline__ void kernel_rot(const size_t i, c64 *a,
	const struct paulis code, const double c, const double s)
{
	c64 sz = s;
	const uint64_t j = paulis_effect(code, i, &sz);
	if (j < i)
		return;

	const c64 xi = a[i], xj = a[j];
	a[i] = c * xi + _Complex_I * conj(sz) * xj;
	a[j] = c * xj + _Complex_I * sz * xi;
}

static __inline__ void kernel_add(
	const size_t i, c64 *restrict a, const c64 *restrict b)
{
	a[i] += b[i];
}

static void qreg_paulirot_hi(struct qreg *reg, struct paulis code_hi, c64 *bm)
{
	paulis_shr(&code_hi, reg->qb_lo);
	const uint64_t rnk_rem = paulis_effect(code_hi, reg->wd.rank, nullptr);

	exch_init(reg, rnk_rem);
	paulis_effect(code_hi, rnk_rem, bm);
	exch_waitall(reg);
}

static void qreg_paulirot_lo(struct qreg *reg, const struct paulis *codes_lo,
	const double *angles, const size_t ncodes, const c64 bm)
{
	for (size_t i = 0; i < reg->namp; i++)
		kernel_mix(i, reg->amp, reg->buf, bm);

	for (size_t k = 0; k < ncodes; k++) {
		const double c = cos(angles[k]), s = sin(angles[k]);
		for (size_t i = 0; i < reg->namp; i++)
			kernel_rot(i, reg->amp, codes_lo[k], c, s);
		for (size_t i = 0; i < reg->namp; i++)
			kernel_rot(i, reg->buf, codes_lo[k], c, -s);
	}

	for (size_t i = 0; i < reg->namp; i++)
		kernel_add(i, reg->amp, reg->buf);
}

void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *phis,
	const size_t ncodes)
{
	c64 bm = 1.0;

	qreg_paulirot_hi(reg, code_hi, &bm);
	qreg_paulirot_lo(reg, codes_lo, phis, ncodes, bm);
}
