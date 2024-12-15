#ifndef QREG_H_INTERNAL
#define QREG_H_INTERNAL

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
