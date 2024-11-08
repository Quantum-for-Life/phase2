#include "c23_compat.h"
#include <complex.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2/paulis.h"

#include <phase2/qreg.h>

struct paulis paulis_new(void)
{
	struct paulis code = {
		.pak = { 0, 0 }
	};

	return code;
}

void paulis_set(struct paulis *code, const int op, const uint32_t n)
{
	const uint64_t n_mask = UINT64_C(1) << n;

	switch (op) {
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
	default:
		break;
	}
}

int paulis_get(const struct paulis code, const uint32_t n)
{
	int pa = (code.pak[0] >> n & 1) | ((code.pak[1] >> n & 1) << 1);

	switch (pa) {
	case 0:
		return PAULI_I;
	case 1:
		return PAULI_X;
	case 2:
		return PAULI_Z; /* sic! */
	case 3:
		return PAULI_Y;
	default:
		unreachable();
	}
}

int paulis_eq(const struct paulis code1, const struct paulis code2)
{
	return code1.pak[0] == code2.pak[0] && code1.pak[1] == code2.pak[1];
}

void paulis_shl(struct paulis *code, const uint32_t n)
{
	code->pak[0] <<= n;
	code->pak[1] <<= n;
}

void paulis_shr(struct paulis *code, const uint32_t n)
{
	code->pak[0] >>= n;
	code->pak[1] >>= n;
}

uint64_t paulis_effect(
	const struct paulis code, const uint64_t i, _Complex double *z)
{
	if (!z)
		goto rt;

	const int mi = stdc_count_ones_ul(i & code.pak[1]);
	const int is = stdc_count_ones_ul(code.pak[0] & code.pak[1]);
	const int r4 = (is + 2 * mi) & 0x3;
	switch (r4) {
	case 0:
		break;
	case 1:
		*z *= I;
		break;
	case 2:
		*z *= -1.0;
		break;
	case 3:
		*z *= -I;
		break;
	default:
		unreachable();
	}

rt:
	return i ^ code.pak[0];
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

void paulis_merge(struct paulis *code, const uint32_t qb_lo,
	const uint32_t qb_hi, const struct paulis lo, const struct paulis hi)
{
	uint64_t mask_lo, mask_hi;

	mask_lo = (UINT64_C(1) << qb_lo) - 1;
	mask_hi = (UINT64_C(1) << (qb_hi + qb_lo)) - 1 - mask_lo;

	code->pak[0] = (lo.pak[0] & mask_lo) | (hi.pak[0] & mask_hi);
	code->pak[1] = (lo.pak[1] & mask_lo) | (hi.pak[1] & mask_hi);
}

int paulis_cmp(const struct paulis a, const struct paulis b)
{
	for (uint32_t n = 0; n < QREG_MAX_WIDTH; n++) {
		const int x = paulis_get(a, n);
		const int y = paulis_get(b, n);
		if (x < y)
			return -1;
		if (x > y)
			return 1;
	}

	return 0;
}