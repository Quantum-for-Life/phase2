#ifndef PAULIS_H
#define PAULIS_H

#include <stdint.h>

typedef _Complex double c64;

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

void paulis_shr(struct paulis *code, uint32_t n);

uint64_t paulis_effect(struct paulis code, uint64_t i, c64 *z);

void paulis_split(struct paulis code, uint32_t qb_lo, uint32_t qb_hi,
	struct paulis *lo, struct paulis *hi);


#endif /* PAULIS_H */
