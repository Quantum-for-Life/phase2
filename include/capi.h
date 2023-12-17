#ifndef PHASE2_CAPI_H
#define PHASE2_CAPI_H

#include <stdlib.h>

#include "circ.h"

void
capi_hadamard(struct circ *c, qbid qb);

void
capi_sgate(struct circ *c, qbid qb);

double
capi_prob0(struct circ *c, qbid qb);

void
capi_blankstate(struct circ *c);

void
capi_set_sysamp(struct circ *c, size_t idx, _Complex double amp);

void
capi_ctl_rotate_pauli(struct circ *c, int *paulis, double angle);

#endif // PHASE2_CAPI_H
