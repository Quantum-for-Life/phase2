#ifndef PAULIS_H
#define PAULIS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if PHASE2_BACKEND == 1 /* QuEST */
#include "QuEST.h"
#else
/* QuEST defines these as enum pauliOpType.
 * To avoid the clash, we pass ints to functions, instead of the enum.
 */
enum {
	PAULI_I = 0,
	PAULI_X = 1,
	PAULI_Y = 2,
	PAULI_Z = 3,
};
#endif /* PHASE2_BACKEND == 1 */

struct paulis {
	uint64_t pak[2];
};

struct paulis paulis_new(void);

int paulis_get(struct paulis code, uint32_t n);

void paulis_set(struct paulis *code, int op, uint32_t n);

int paulis_eq(struct paulis code1, struct paulis code2);

void paulis_shl(struct paulis *code, uint32_t n);

void paulis_shr(struct paulis *code, uint32_t n);

uint64_t paulis_effect(struct paulis code, uint64_t i, _Complex double *z);

void paulis_split(struct paulis code, uint32_t qb_lo, uint32_t qb_hi,
	struct paulis *lo, struct paulis *hi);

void paulis_merge(struct paulis *code, uint32_t qb_lo, uint32_t qb_hi,
	struct paulis lo, struct paulis hi);

#ifdef __cplusplus
}
#endif

#endif /* PAULIS_H */
