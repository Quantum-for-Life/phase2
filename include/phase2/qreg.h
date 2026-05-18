#ifndef QREG_H
#define QREG_H

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

int qreg_init(struct qreg *reg, uint32_t qb);

void qreg_free(struct qreg *reg);

void qreg_getamp(struct qreg *reg, uint64_t i, _Complex double *z);

void qreg_setamp(struct qreg *reg, uint64_t i, _Complex double z);

void qreg_zero(struct qreg *reg);

/* Apply product of Pauli rotations for a list of strings sharing the same
 * operation on hi qubits.
 */
void qreg_paulirot(struct qreg *reg, struct paulis code_hi,
	const struct paulis *codes_lo, const double *phis, size_t ncodes);

#ifdef __cplusplus
}
#endif

#endif // QREG_H
