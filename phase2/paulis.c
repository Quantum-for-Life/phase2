#include <complex.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2/paulis.h"

typedef _Complex double c64;

struct paulis paulis_new(void)
{
	struct paulis code = {
		.pak = { 0, 0 }
	};

	return code;
}

static int paulis_countis(struct paulis code)
{
	return __builtin_popcountll(code.pak[0] & code.pak[1]);
}

void paulis_set(
	struct paulis *code, const enum pauli_op pauli, const uint32_t n)
{
	const uint64_t n_mask = UINT64_C(1) << n;

	switch (pauli) {
	case PAULI_I:
		code->pak[0] &= ~n_mask;
		code->pak[1] &= ~n_mask;
		break;
	case PAULI_X:
		code->pak[0] |= n_mask;
		code->pak[1] &= ~n_mask;
		break;
	case PAULI_Y:
		code->pak[0] |= n_mask;
		code->pak[1] |= n_mask;
		break;
	case PAULI_Z:
		code->pak[0] &= ~n_mask;
		code->pak[1] |= n_mask;
		break;
	}
}

enum pauli_op paulis_get(const struct paulis code, const uint32_t n)
{
	int pa = 0;
	pa |= code.pak[0] >> n & 1;
	pa |= (code.pak[1] >> n & 1) << 1;

	switch (pa) {
	case 0:
		return PAULI_I;
	case 1:
		return PAULI_X;
	case 2:
		return PAULI_Z;
	case 3:
		return PAULI_Y;
	default:
		__builtin_unreachable();
	}
}

int paulis_eq(const struct paulis code1, const struct paulis code2)
{
	return code1.pak[0] == code2.pak[0] && code1.pak[1] == code2.pak[1];
}

void paulis_shr(struct paulis *code, const uint32_t n)
{
	code->pak[0] >>= n;
	code->pak[1] >>= n;
}

uint64_t paulis_effect(const struct paulis code, const uint64_t i, c64 *z)
{
	uint64_t j = i ^ code.pak[0];
	if (z != NULL) {
		const int minus = __builtin_popcountll(j & code.pak[1]);
		const int root4 = (paulis_countis(code) + 2 * minus) & 0x3;

		switch (root4) {
		case 0:
			break;
		case 1:
			*z *= -I;
			break;
		case 2:
			*z *= -1.0;
			break;
		case 3:
			*z *= I;
			break;
		default:
			__builtin_unreachable();
		}
	}

	return j;
}

void paulis_split(const struct paulis code, const uint32_t qb_lo,
	const uint32_t qb_hi, struct paulis *lo, struct paulis *hi)
{
	uint64_t mask_lo, mask_hi;

	mask_lo = (UINT64_C(1) << qb_lo) - 1;
	mask_hi = (UINT64_C(1) << (qb_hi + qb_lo)) - 1 - mask_lo;

	lo->pak[0] = code.pak[0] & mask_lo;
	lo->pak[1] = code.pak[1] & mask_lo;

	hi->pak[0] = code.pak[0] & mask_hi;
	hi->pak[1] = code.pak[1] & mask_hi;
}
