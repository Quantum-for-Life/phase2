#ifndef PAULIS_H
#define PAULIS_H

#include "types.h"

#define PAULI_MAX_WIDTH (64)

typedef enum root4 {
	R0, // +1
	R1, // +i
	R2, // -1
	R3, // -i:
} root4;

enum pauli_op {
	PAULI_I = 0,
	PAULI_X = 1,
	PAULI_Y = 2,
	PAULI_Z = 3,
};

static const char PAULI_LABEL[4] = { 'I', 'X', 'Y', 'Z' };

struct paulis {
	u64 pak[2];
};

struct paulis paulis_new(void);
struct paulis paulis_from_ops(const enum pauli_op *paulis, u32 num_paulis);

void paulis_print(struct paulis code);

void	      paulis_set(struct paulis *code, enum pauli_op pauli, u32 n);
enum pauli_op paulis_get(struct paulis code, u32 n);

int paulis_eq(struct paulis code1, struct paulis code2);

void paulis_mask(struct paulis *code, u64 mask);
void paulis_shr(struct paulis *code, u32 n);
u64  paulis_effect(struct paulis code, u64 i, root4 *z);
void paulis_split(struct paulis code, u32 qb_lo, u32 qb_hi, struct paulis *lo,
	struct paulis *hi);

#endif // PAULIS_H
