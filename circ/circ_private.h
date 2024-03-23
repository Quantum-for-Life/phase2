#ifndef PHASE2_CIRC_PRIVATE_H
#define PHASE2_CIRC_PRIVATE_H

#include "qreg.h"

#include "circ.h"

struct circ {
	struct circuit *ct;
	void	     *data;
	int	    *cl, *qb;

	/* Qubit register */
	struct qreg quest_qureg;
};

#endif // PHASE2_CIRC_PRIVATE_H
