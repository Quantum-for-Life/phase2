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
	reg->backend = nullptr;

	return 0;
}

void qreg_backend_free(struct qreg *reg)
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
	MPI_Barrier(MPI_COMM_WORLD);
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

/*
 * kernel_mix - project local and partner amplitudes into the
 * +/-1 eigenspaces of the hi-part Pauli operator.
 *
 * After MPI exchange, a[] holds this rank's amplitudes and
 * b[] holds the partner rank's amplitudes.  The hi-part
 * Pauli acting on the partner rank's index produces phase
 * bm (a fourth root of unity).
 *
 * The mixing step computes:
 *   a[i] = (a[i] + bm*b[i]) / 2
 *   b[i] = (a[i] - bm*b[i]) / 2
 *
 * This is the projector decomposition:
 *   (I + bm*P_hi)/2  projects into the +1 eigenspace,
 *   (I - bm*P_hi)/2  projects into the -1 eigenspace.
 *
 * After mixing, rotations are applied independently to
 * each half, then the halves are recombined via addition.
 */
static __inline__ void kernel_mix(
	const size_t i, c64 *restrict a, c64 *restrict b, const c64 bm)
{
	b[i] *= bm;

	const c64 x = a[i];
	const c64 y = b[i];
	a[i] = (x + y) / 2.0;
	b[i] = (x - y) / 2.0;
}

/*
 * kernel_rot - apply Pauli rotation exp(i*phi*P_lo) to the
 * local amplitude vector.
 *
 * For each basis index i, compute P_lo|i> = z|j>.  The
 * rotation couples the pair (i, j) as a 2x2 block:
 *   a[i] <- cos(phi)*a[i] + i*conj(z)*sin(phi)*a[j]
 *   a[j] <- cos(phi)*a[j] + i*z*sin(phi)*a[i]
 *
 * The guard j < i ensures each (i, j) pair is processed
 * exactly once: when i is the smaller index of the pair.
 * If P_lo|i> = z|i> (diagonal case, j == i), the guard
 * still works since j < i is false and the rotation
 * reduces to a phase.
 */
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

/*
 * qreg_paulirot_hi - MPI amplitude exchange for the hi-part
 * Pauli operator.
 *
 * The hi qubits are distributed across MPI ranks.  Applying
 * paulis_effect to this rank's index gives the partner rank
 * that holds the paired amplitudes.  Non-blocking send/recv
 * exchanges the full local amplitude vector with the partner.
 * The phase bm from applying P_hi to the partner's rank
 * index is returned for use in kernel_mix.
 */
static void qreg_paulirot_hi(struct qreg *reg, struct paulis code_hi, c64 *bm)
{
	paulis_shr(&code_hi, reg->qb_lo);
	const uint64_t rnk_rem = paulis_effect(code_hi, reg->wd.rank, nullptr);

	exch_init(reg, rnk_rem);
	paulis_effect(code_hi, rnk_rem, bm);
	exch_waitall(reg);
}

/*
 * qreg_paulirot_lo - apply batched lo-qubit Pauli rotations.
 *
 * Three-phase protocol:
 *  1. Mix: project amp[] and buf[] (partner data) into the
 *     +/-1 eigenspaces of the hi-part operator using bm.
 *  2. Rotate: apply each lo-part rotation to both halves.
 *     amp[] gets +phi, buf[] gets -phi (opposite sign),
 *     so the eigenspaces evolve independently.
 *  3. Add: recombine amp[] += buf[] to reconstruct the
 *     full state vector on this rank.
 */
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
	const struct paulis *codes_lo, const double *phis, const size_t ncodes)
{
	c64 bm = 1.0;

	qreg_paulirot_hi(reg, code_hi, &bm);
	qreg_paulirot_lo(reg, codes_lo, phis, ncodes, bm);
}
