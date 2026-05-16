#ifndef PAULIS_H
#define PAULIS_H

/*
 * phase2/paulis.h -- bit-packed Pauli-string
 * encoding and the operations on it.
 *
 * A Pauli string acts on up to QREG_MAX_WIDTH = 64
 * qubits, with one Pauli operator (I, X, Y, or Z)
 * per qubit.  Each string is stored in two 64-bit
 * words using the standard symplectic
 * representation (see paulis.c for the derivation):
 *
 *     pak[0]  =  flip bits   (X / Y components)
 *     pak[1]  =  phase bits  (Z / Y components)
 *
 * Qubit n is at bit n of each word.  The per-qubit
 * code is read off via paulis_get / set:
 *
 *     I = 00,  X = 10,  Z = 01,  Y = 11
 *
 * All public operations take qubit indices in
 * [0, 64); higher indices are undefined (the shift
 * `<< n` would wrap).  The struct itself is 16
 * bytes and cheap to pass by value -- nearly every
 * function does.
 */

#include <stdint.h>

/*
 * Single-qubit Pauli operators.  Integer values
 * are stable and form part of the public API;
 * paulis_set encodes them into the symplectic form
 * via a switch, paulis_get inverts.
 */
enum pauli_op {
	PAULI_I = 0,
	PAULI_X = 1,
	PAULI_Y = 2,
	PAULI_Z = 3,
};

/*
 * Bit-packed n-qubit Pauli string.
 *
 *   pak[0]   flip / X-component bits; bit n encodes
 *            qubit n.
 *   pak[1]   phase / Z-component bits; bit n encodes
 *            qubit n.
 *
 * (pak[0]:n, pak[1]:n) read off as (I, X, Z, Y) for
 * (00, 10, 01, 11) respectively.  Default-zero
 * initialiser is the all-identity string.
 */
struct paulis {
	uint64_t pak[2];
};

/* Return the all-identity Pauli string. */
struct paulis paulis_new(void);

/*
 * Decode the Pauli operator at qubit n
 * (0 <= n < 64).  Returns one of PAULI_{I,X,Y,Z}.
 */
enum pauli_op paulis_get(struct paulis code, uint32_t n);

/*
 * Encode the Pauli operator at qubit n
 * (0 <= n < 64).  Replaces whatever was there
 * before; bits at other qubits are preserved.
 * `op` must be one of PAULI_{I,X,Y,Z} -- any other
 * value leaves the string unchanged.
 */
void paulis_set(struct paulis *code, enum pauli_op op, uint32_t n);

/*
 * Bitwise equality of two Pauli strings.  Returns
 * 1 if every qubit's Pauli is the same, 0
 * otherwise.
 */
int paulis_eq(struct paulis code1, struct paulis code2);

/*
 * Shift the encoded Pauli string n positions
 * towards higher / lower qubit indices.  Bits
 * shifted out are dropped (undefined-behaviour
 * free for n in [0, 64)).  No bounds checking on
 * the resulting effective qubit count -- callers
 * that mix shift and split must track the active
 * range themselves.
 */
void paulis_shl(struct paulis *code, uint32_t n);
void paulis_shr(struct paulis *code, uint32_t n);

/*
 * Apply the Pauli string to computational basis
 * state |i>: P|i> = z * |j>, returning j and
 * multiplying *z by the phase factor.
 *
 * `z` may be NULL, in which case the phase is not
 * computed and only the flipped index j is
 * returned.  The cold path: callers that need only
 * the partner index (e.g. qreg_paulirot_hi
 * computing the partner MPI rank) pass NULL to
 * skip two popcounts and a branch.
 */
uint64_t paulis_effect(struct paulis code, uint64_t i, _Complex double *z);

/*
 * Integer core of paulis_effect, shared between the
 * host implementation (phase2/paulis.c) and the
 * CUDA device kernel (phase2/qreg_cuda_lo.cu).
 *
 * Returns j = i XOR pak[0] (the flipped basis-state
 * index) and writes the phase exponent r4 in
 * [0, 4) to *r4_out: the amplitude on |j> picks up
 * a factor of (i)^r4_out, where i is the imaginary
 * unit.  The caller does the complex multiply in
 * its own complex type (host uses _Complex double,
 * device uses cuDoubleComplex).
 *
 * Any change to the phase formula belongs here and
 * here only -- the host and device wrappers both
 * read r4_out and translate.
 *
 * `__host__ __device__` under nvcc; plain inline
 * otherwise.
 */
#ifdef __CUDACC__
#  define PAULIS_HD __host__ __device__
#else
#  define PAULIS_HD
#endif

PAULIS_HD static inline uint64_t paulis_effect_raw(struct paulis code,
	uint64_t i, int *r4_out)
{
	const int mi = __builtin_popcountll(i & code.pak[1]);
	const int is = __builtin_popcountll(code.pak[0] & code.pak[1]);
	*r4_out = (is + 2 * mi) & 0x3;
	return i ^ code.pak[0];
}

/*
 * Split `code` into two disjoint qubit ranges:
 *
 *   *lo  carries qubits [0, qb_lo).
 *   *hi  carries qubits [qb_lo, qb_lo + qb_hi).
 *
 * Bits outside those two ranges are zeroed in the
 * outputs.  Used by the rotation cache to separate
 * MPI-distributed hi qubits from local lo qubits.
 * Inputs and outputs may not alias.
 */
void paulis_split(struct paulis code, uint32_t qb_lo, uint32_t qb_hi,
	struct paulis *lo, struct paulis *hi);

/*
 * Lexicographic total order on Pauli strings,
 * highest qubit first.  Returns -1, 0, or +1.  Used
 * by circ_hamil_sort_lex to group Hamiltonian
 * terms with matching hi-qubit codes contiguously,
 * which maximises the batch-cache hit rate.
 */
int paulis_cmp(struct paulis a, struct paulis b);

#endif /* PAULIS_H */
