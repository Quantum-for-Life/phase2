/**
 * Circuit interface.
 */
#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

struct circ {
	struct circuit *ct;
	void *data;

	/* Qubit register */
	void *reg;

	int *mea_cl;
	int *mea_qb;
	int *sys_qb;
	int *anc_qb;
};

typedef int (*circ_op)(struct circ *);

struct circuit {
	const char *name;
	void *data;

	size_t num_mea_qb;
	size_t num_sys_qb;
	size_t num_anc_qb;

	circ_op reset;
	circ_op prepst;
	circ_op effect;
	circ_op measure;
};

int circ_init(struct circ *c, struct circuit *ct, void *data);

void circ_destroy(struct circ *c);

void circ_report(struct circ const *c);

int circ_reset(struct circ *c);

int circ_simulate(struct circ *c);

#endif //PHASE2_CIRC_H
