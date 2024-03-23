#ifndef QREG_H
#define QREG_H

#include <mpi.h>

#include "paulis.h"
#include "types.h"

#define QREG_MAX_NUM_QUBITS (64)

struct qreg {
	const struct ev *ev;

	u32 qb_lo, qb_hi;

	fl *amp[2], *buf[2];
	u64  num_amps;

	int	     msg_count;
	MPI_Request *reqs_snd, *reqs_rcv;
	size_t	     num_reqs;
};

int  qreg_init(struct qreg *reg, u32 num_qubits, const struct ev *ev);
void qreg_destroy(struct qreg *reg);

void qreg_getamp(const struct qreg *reg, u64 n, fl (*z)[2]);
void qreg_setamp(struct qreg *reg, u64 n, const fl z[2]);

void qreg_zerostate(struct qreg *reg);
// void qreg_paulirot(struct qreg *reg, struct paulis code, f64 angle);
void qreg_paulirot(struct qreg *reg, struct paulis code_hi,
	const struct paulis *codes_lo, const fl *angles, size_t num_codes);

#endif // QREG_H
