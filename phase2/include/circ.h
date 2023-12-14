/**
 * Circuit interface.
 */
#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

typedef size_t qbid;

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

struct circ *circ_init(struct circuit *ct, void *data);

void circ_destroy(struct circ *c);

void circ_report(struct circ const *c);

int circ_reset(struct circ *c);

int circ_simulate(struct circ *c);

qbid circ_mea_qb(struct circ *c, size_t idx);

qbid circ_sys_qb(struct circ *c, size_t idx);

qbid circ_anc_qb(struct circ *c, size_t idx);

/* Quantum API */
void circ_hadamard(struct circ *c, qbid q);

void circ_sgate(struct circ *c, qbid q);

double circ_prob0(struct circ *c, qbid q);

void circ_blankstate(struct circ *c);

void circ_setsysamp(struct circ *c, size_t idx, _Complex double amp);

void circ_sys_control_rotate_pauli(struct circ *c, int *paulis, double angle);

#endif //PHASE2_CIRC_H
