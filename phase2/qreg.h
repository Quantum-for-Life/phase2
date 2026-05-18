#ifndef QREG_H_INTERNAL
#define QREG_H_INTERNAL

/*
 * Subsystem-private surface shared between the
 * dispatcher (qreg.c) and the backends (qreg_qreg.c
 * CPU, qreg_cuda.c GPU).  The paulirot dispatch and
 * the file-local exchange helpers live duplicated
 * in each backend deliberately: hoisting them into
 * qreg.c crosses a translation-unit boundary that
 * defeats inlining and adds a measurable per-call
 * regression on b-qreg / b-algos.  Touch with
 * `make bench` running.
 */

#ifdef __cplusplus
extern "C" {
#endif

uint64_t qreg_getilo(const struct qreg *reg, uint64_t i);

uint64_t qreg_getihi(const struct qreg *reg, uint64_t i);

int qreg_backend_init(struct qreg *reg);

void qreg_backend_free(struct qreg *reg);

#ifdef __cplusplus
}
#endif

#endif /* QREG_H_INTERNAL */
