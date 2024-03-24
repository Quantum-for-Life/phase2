#include <stdio.h>

#include "paulis.h"

struct paulis paulis_new(void)
{
	struct paulis code;
	code.pak[0] = (u64)0;
	code.pak[1] = (u64)0;

	return code;
}

struct paulis paulis_from_ops(const enum pauli_op *paulis, const u32 num_paulis)
{
	struct paulis code = paulis_new();
	for (u32 i = 0; i < num_paulis; i++)
		paulis_set(&code, paulis[i], i);

	return code;
}

void paulis_print(const struct paulis code)
{
	for (int i = 0; i < PAULI_MAX_WIDTH; i++) {
		const enum pauli_op pauli = paulis_get(code, i);
		printf("%c", PAULI_LABEL[pauli]);
	}
	printf("\n");
}

void paulis_set(struct paulis *code, const enum pauli_op pauli, const u32 n)
{
	const u64 n_mask = (u64)1 << n;

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

enum pauli_op paulis_get(const struct paulis code, const u32 n)
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

void paulis_mask(struct paulis *code, const u64 mask)
{
	code->pak[0] &= mask;
	code->pak[1] &= mask;
}

void paulis_shr(struct paulis *code, const u32 n)
{
	code->pak[0] >>= n;
	code->pak[1] >>= n;
}

u64 paulis_effect(struct paulis code, u64 i, root4 *z)
{
	const u64 *p = code.pak;
	const u64  j = i ^ p[0];

	if (z == NULL)
		return j;

	int root4 = 0;
	for (u64 yp = ~i & p[0] & p[1]; yp > 0; yp >>= 1)
		if (yp & 1)
			root4 += 1;
	for (u64 zm = i & ~p[0] & p[1]; zm > 0; zm >>= 1)
		if (zm & 1)
			root4 += 2;
	for (u64 ym = i & p[0] & p[1]; ym > 0; ym >>= 1)
		if (ym & 1)
			root4 += 3;
	*z = root4 & 3;

	return j;
}

void paulis_split(const struct paulis code, const u32 qb_lo, const u32 qb_hi,
	struct paulis *lo, struct paulis *hi)
{
	const u64 mask_lo = ((u64)1 << qb_lo) - 1;

	*lo = code;
	paulis_mask(lo, mask_lo);

	*hi = code;
	paulis_mask(hi, ((u64)1 << (qb_lo + qb_hi)) - mask_lo);
}
