#ifndef QREG_H
#define QREG_H

#include <stdint.h>

#include "mpi.h"

#include "paulis.h"

#define QREG_MAX_WIDTH (64)

struct qreg {
	struct qreg_ev {
		int num_ranks;
		int rank;
	} ev;

	uint32_t qb_lo, qb_hi;

	c64 *amp, *buf;
	uint64_t num_amps;

	int msg_count;
	MPI_Request *reqs_snd, *reqs_rcv;
	size_t num_reqs;
};

int qreg_init(struct qreg *reg, uint32_t num_qubits);

void qreg_destroy(struct qreg *reg);

void qreg_getamp(const struct qreg *reg, uint64_t i, c64 *z);

void qreg_setamp(struct qreg *reg, uint64_t i, c64 z);

void qreg_zero(struct qreg *reg);

void qreg_paulirot(struct qreg *reg, struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles, size_t num_codes);

#endif // QREG_H
