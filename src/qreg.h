#ifndef QREG_H
#define QREG_H

#include <stdint.h>

#include <mpi.h>

#define PAULI_MAX_WIDTH (64)

typedef int8_t	i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef uint8_t	 u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float  f32;
typedef double f64;

typedef f64 fl;

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
enum pauli_op paulis_get(struct paulis code, u32 n);

void paulis_set(struct paulis *code, enum pauli_op pauli, u32 n);
int  paulis_eq(struct paulis code1, struct paulis code2);
void paulis_mask(struct paulis *code, u64 mask);
void paulis_shr(struct paulis *code, u32 n);
u64  paulis_effect(struct paulis code, u64 i, root4 *z);
void paulis_split(struct paulis code, u32 qb_lo, u32 qb_hi, struct paulis *lo,
	struct paulis *hi);

struct qreg {
	struct qreg_ev {
		int num_ranks;
		int rank;
	} ev;

	u32 qb_lo, qb_hi;

	fl *amp[2], *buf[2];
	u64 num_amps;

	int	     msg_count;
	MPI_Request *reqs_snd, *reqs_rcv;
	size_t	     num_reqs;
};

int  qreg_init(struct qreg *reg, u32 num_qubits);
void qreg_destroy(struct qreg *reg);

void qreg_getamp(const struct qreg *reg, u64 n, fl (*z)[2]);
void qreg_setamp(struct qreg *reg, u64 n, const fl z[2]);

void qreg_zero(struct qreg *reg);
void qreg_paulirot(struct qreg *reg, struct paulis code_hi,
	const struct paulis *codes_lo, const fl *angles, size_t num_codes);

#endif // QREG_H
