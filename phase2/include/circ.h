/**
 * Circuit interface.
 */
#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

/* Dynamical instance of a circuit */
struct circ;

typedef size_t qbid;

struct circuit {
	const char *name;

	size_t num_mea_qb;
	size_t num_sys_qb;
	size_t num_anc_qb;

	int (*reset)(struct circ *);
	int (*prepst)(struct circ *);
	int (*effect)(struct circ *);
	int (*measure)(struct circ *);
};

int circ_env_initialize();

int circ_env_shutdown();

struct circ *circ_init(struct circuit *ct, void *data);

void circ_destroy(struct circ *c);

void *circ_data(struct circ *c);

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
