#ifndef PHASE2_CIRC_INTL_H
#define PHASE2_CIRC_INTL_H

#include "QuEST.h"

#include "circ.h"

Qureg
circ_intl_quest_qureg(const struct circ *c);

int *
circ_intl_get_qb(const struct circ *c);

#endif // PHASE2_CIRC_INTL_H
