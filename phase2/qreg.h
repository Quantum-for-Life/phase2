#ifndef QREG_H_INTERNAL
#define QREG_H_INTERNAL

/*
 * Subsystem-private surface shared between the
 * dispatcher (qreg.c) and the backends (qreg_qreg.c
 * CPU, qreg_cuda.c GPU).
 */

#ifdef __cplusplus
extern "C" {
#endif

uint64_t qreg_getilo(const struct qreg *reg, uint64_t i);

uint64_t qreg_getihi(const struct qreg *reg, uint64_t i);

int qreg_backend_init(struct qreg *reg);

void qreg_backend_free(struct qreg *reg);

/* Post one MPI exchange round: reg's local amp -> rnk_rem,
 * rnk_rem's local amp -> reg's buf. */
void qreg_backend_exch_init(struct qreg *reg, int rnk_rem);

void qreg_backend_exch_waitall(struct qreg *reg);

/* Apply ncodes lo-qubit Pauli rotations to the post-
 * exchange buffers; bm is the hi-Pauli phase carried
 * from the partner rank. */
void qreg_backend_paulirot_lo(struct qreg *reg,
	const struct paulis *codes_lo, const double *angles,
	size_t ncodes, _Complex double bm);

#ifdef __cplusplus
}
#endif

#endif /* QREG_H_INTERNAL */
