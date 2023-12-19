#ifndef PHASE2_CIRC_PRIVATE_H
#define PHASE2_CIRC_PRIVATE_H

#include "QuEST.h"

#include "circ.h"

struct circ {
	struct circuit *ct;
	void	       *data;
	int	       *cl, *qb;

	/* Qubit register */
	Qureg quest_qureg;
};

#endif // PHASE2_CIRC_PRIVATE_H
