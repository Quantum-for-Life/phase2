#ifndef QREG_H
#define QREG_H

#include <stdint.h>

#include "mpi.h"

#define QREG_MAX_WIDTH (64)

enum pauli_op {
	PAULI_I = 0,
	PAULI_X = 1,
	PAULI_Y = 2,
	PAULI_Z = 3,
};

static const char PAULI_LABEL[4] = { 'I', 'X', 'Y', 'Z' };

struct paulis {
	uint64_t pak[2];
};

struct paulis paulis_new(void);

enum pauli_op paulis_get(struct paulis code, uint32_t n);

void paulis_set(struct paulis *code, enum pauli_op pauli, uint32_t n);

int paulis_eq(struct paulis code1, struct paulis code2);

void paulis_mask(struct paulis *code, uint64_t mask);

void paulis_shr(struct paulis *code, uint32_t n);

uint64_t paulis_effect(struct paulis code, uint64_t i, _Complex double *z);

void paulis_split(struct paulis code, 
	uint32_t qb_lo, uint32_t qb_hi, struct paulis *lo,struct paulis *hi);

struct qreg {
	struct qreg_ev {
		int num_ranks;
		int rank;
	} ev;

	uint32_t qb_lo, qb_hi;

	_Complex double *amp;
	_Complex double *buf;
	uint64_t num_amps;

	int	     msg_count;
	MPI_Request *reqs_snd, *reqs_rcv;
	size_t	     num_reqs;
};

int qreg_init(struct qreg *reg, uint32_t num_qubits);

void qreg_destroy(struct qreg *reg);

void qreg_getamp(const struct qreg *reg, uint64_t i, _Complex double *z);

void qreg_setamp(struct qreg *reg, uint64_t i, _Complex double z);

void qreg_zero(struct qreg *reg);

void qreg_paulirot(struct qreg *reg, struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles, size_t num_codes);

#endif // QREG_H
