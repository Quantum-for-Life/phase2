#ifndef CIRC_IMPL_H
#define CIRC_IMPL_H

#include "phase2/circ.h"

#ifdef __cplusplus
extern "C" {
#endif

int circ_res_init(struct circ *c);

void circ_res_destroy(struct circ *c);

#ifdef __cplusplus
}
#endif

#endif // CIRC_IMPL_H
