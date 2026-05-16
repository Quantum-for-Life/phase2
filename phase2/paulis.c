#include "c23_compat.h"
#include <complex.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"

inline struct paulis paulis_new(void)
{
	struct paulis code = {
		.pak = { 0, 0 }
	};

	return code;
}

/*
 * paulis_set - encode a single-qubit Pauli operator into the
 * bit-packed representation.
 *
 * Each n-qubit Pauli string is stored as two 64-bit words:
 *   pak[0] = "flip bits"  (X component)
 *   pak[1] = "phase bits" (Z component)
 *
 * Bit n of each word encodes qubit n.  The mapping is:
 *   I = 00   (pak[0]:0, pak[1]:0)
 *   X = 10   (pak[0]:1, pak[1]:0)
 *   Z = 01   (pak[0]:0, pak[1]:1)
 *   Y = 11   (pak[0]:1, pak[1]:1)   (Y = iXZ)
 *
 * This is the standard "symplectic" representation used in
 * stabiliser / Pauli-group literature.
 */
void paulis_set(struct paulis *code, enum pauli_op op, uint32_t n)
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

/*
 * paulis_get - decode a single-qubit Pauli from the packed form.
 *
 * The two-bit value pa = (pak[1][n] << 1) | pak[0][n] yields:
 *   0 -> I,  1 -> X,  2 -> Z,  3 -> Y
 *
 * Note: case 2 (binary 10) maps to Z, not Y.  The internal
 * bit layout (flip, phase) means pak[1]=1, pak[0]=0 is Z,
 * which reads as the integer 2 — easily mistaken for Y if
 * one assumes a naive I/X/Y/Z = 0/1/2/3 encoding.
 */
enum pauli_op paulis_get(struct paulis code, uint32_t n)
{
	const int pa = (code.pak[0] >> n & 1) | ((code.pak[1] >> n & 1) << 1);

	switch (pa) {
	case 0:
		return PAULI_I;
	case 1:
		return PAULI_X;
	case 2:
		return PAULI_Z; /* see comment above */
	case 3:
		return PAULI_Y;
	default:
		unreachable();
	}
}

inline int paulis_eq(struct paulis code1, struct paulis code2)
{
	return code1.pak[0] == code2.pak[0] && code1.pak[1] == code2.pak[1];
}

inline void paulis_shl(struct paulis *code, uint32_t n)
{
	code->pak[0] <<= n;
	code->pak[1] <<= n;
}

inline void paulis_shr(struct paulis *code, uint32_t n)
{
	code->pak[0] >>= n;
	code->pak[1] >>= n;
}

/*
 * paulis_effect - compute the action of a Pauli string on a
 * computational basis state:  P|i> = z * |j>.
 *
 * Output:
 *   returns j = i XOR pak[0]   (X/Y bits flip the state)
 *   multiplies *z by the phase factor i^r, where
 *     r = (is + 2*mi) mod 4
 *
 * Derivation:
 *   Write P = (product of single-qubit Paulis).  Each Y_k
 *   contributes a factor of i (since Y = iXZ), and each Z_k
 *   or Y_k acting on a |1> contributes a factor of -1.
 *
 *   is = popcount(pak[0] & pak[1])
 *      = number of Y factors, each giving i^1.
 *   mi = popcount(i & pak[1])
 *      = number of qubits where the state bit is 1 AND the
 *        operator has a Z or Y component, each giving (-1).
 *
 *   Total phase = i^is * (-1)^mi = i^(is + 2*mi).
 *   Reduce mod 4 to index into {1, i, -1, -i}.
 *
 * Thin wrapper around paulis_effect_raw (paulis.h),
 * which is the integer core shared with the CUDA
 * device kernel.  Only the _Complex double multiply
 * lives here.
 */
uint64_t paulis_effect(struct paulis code, uint64_t i, _Complex double *z)
{
	int r4;
	const uint64_t j = paulis_effect_raw(code, i, &r4);
	if (!z)
		return j;

	/* Switch (not a table indexed by r4) so the case-0
	 * path skips the complex multiply entirely.  The
	 * device kernel uses the table form: GPU warps
	 * prefer the uniform load over a switch that would
	 * diverge across threads. */
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
	return j;
}

void paulis_split(struct paulis code, uint32_t qb_lo, uint32_t qb_hi,
	struct paulis *lo, struct paulis *hi)
{
	/* mask_lo covers bits [0, qb_lo); mask_hi
	 * covers bits [qb_lo, qb_lo + qb_hi).  Bits
	 * outside that combined range fall through to
	 * neither output. */
	const uint64_t mask_lo = (UINT64_C(1) << qb_lo) - 1;
	const uint64_t mask_hi = (UINT64_C(1) << (qb_hi + qb_lo)) - 1 - mask_lo;

	lo->pak[0] = code.pak[0] & mask_lo;
	lo->pak[1] = code.pak[1] & mask_lo;

	hi->pak[0] = code.pak[0] & mask_hi;
	hi->pak[1] = code.pak[1] & mask_hi;
}

/*
 * paulis_cmp - lexicographic order, highest qubit
 * first.  The early-out via paulis_eq lets the
 * common "same hi-code" check in the batch cache
 * skip the per-qubit loop.  The loop reads qubits
 * MSB-to-LSB so the resulting order matches the
 * bit-string read top-down.
 */
int paulis_cmp(struct paulis a, struct paulis b)
{
	if (paulis_eq(a, b))
		return 0;

	for (uint32_t n = 0; n < QREG_MAX_WIDTH; n++) {
		const enum pauli_op x = paulis_get(a, QREG_MAX_WIDTH - n - 1);
		const enum pauli_op y = paulis_get(b, QREG_MAX_WIDTH - n - 1);
		if (x < y)
			return -1;
		if (x > y)
			return 1;
	}

	unreachable();
}
