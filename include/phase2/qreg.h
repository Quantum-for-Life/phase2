#ifndef QREG_H
#define QREG_H

#include <stdint.h>

#include "mpi.h"

#include "phase2/paulis.h"
#include "phase2/world.h"

#define QREG_MAX_WIDTH (64)

struct qreg {
	uint32_t qb_lo, qb_hi;

	_Complex double *amp, *buf;
	uint64_t num_amps;

	int msg_count;
	MPI_Request *reqs_snd, *reqs_rcv;
	size_t num_reqs;
};

int qreg_init(struct qreg *reg, uint32_t num_qubits);

void qreg_destroy(struct qreg *reg);

void qreg_getamp(const struct qreg *reg, uint64_t i, _Complex double *z);

void qreg_setamp(struct qreg *reg, uint64_t i, _Complex double z);

void qreg_zero(struct qreg *reg);

/* Apply product of Pauli rotations for a list of strings sharing the same 
 * operation on hi qubits.
 *
 */
void qreg_paulirot(struct qreg *reg, struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles, size_t num_codes);

#endif // QREG_H
