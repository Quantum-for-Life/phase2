#ifndef QREG_H
#define QREG_H

/*
 * MPI-distributed quantum register.  Total qubits =
 * qb_hi + qb_lo, where qb_hi = log2(world size) and
 * each rank owns a slice of 2^qb_lo amplitudes.
 * Operations on the hi qubits cross ranks via MPI;
 * lo-qubit operations stay local.
 */

#include <stdint.h>

#include "mpi.h"

#include "phase2/paulis.h"
#include "phase2/world.h"

#ifdef __cplusplus
extern "C" {
#endif

#define QREG_MAX_WIDTH (64)

struct qreg {
	struct world_info wd;

	uint32_t qb_lo, qb_hi;

	_Complex double *amp, *buf;
	uint64_t namp;

	int msg_count;
	MPI_Request *reqs_snd, *reqs_rcv;
	size_t nreqs;

	/* Backend-private handle.  Opaque to the
	 * public surface; set by qreg_backend_init,
	 * dereferenced inside the backend-specific
	 * sources (qreg_qreg.c, qreg_cuda.c, etc.). */
	void *backend;
};

/* Allocate the rank-local amplitude buffer and the
 * partner-exchange scratch.  Requires
 * qb > log2(world size) and qb <= QREG_MAX_WIDTH;
 * returns -1 otherwise. */
int qreg_init(struct qreg *reg, uint32_t qb);

void qreg_free(struct qreg *reg);

/* Read amplitude at global index i.  Bcast from the
 * owning rank into z on every rank. */
void qreg_getamp(struct qreg *reg, uint64_t i, _Complex double *z);

/* Write amplitude at global index i on the owning
 * rank; all ranks synchronise via barrier. */
void qreg_setamp(struct qreg *reg, uint64_t i, _Complex double z);

/* Zero every amplitude on every rank. */
void qreg_zero(struct qreg *reg);

/* Apply product of Pauli rotations sharing one hi-
 * qubit code: exp(i * phis[k] * code_hi (x) codes_lo[k])
 * for k = 0..ncodes-1.  One MPI exchange amortises
 * across all ncodes rotations. */
void qreg_paulirot(struct qreg *reg, struct paulis code_hi,
	const struct paulis *codes_lo, const double *phis, size_t ncodes);

#ifdef __cplusplus
}
#endif

#endif // QREG_H
