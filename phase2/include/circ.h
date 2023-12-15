/**
 * Circuit interface.
 */
#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

typedef size_t qbid;

struct circ;

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

int circ_initialize();

int circ_shutdown();

struct circ *circ_create(struct circuit *ct, void *data);

void circ_destroy(struct circ *c);

void *circ_data(const struct circ *c);

int circ_report(struct circ const *c);

int circ_reset(struct circ *c);

int circ_run(struct circ *c);

size_t circ_num_meaqb(const struct circ *c);

size_t circ_num_sysqb(const struct circ *c);

size_t circ_num_ancqb(const struct circ *c);

qbid circ_meaqb(const struct circ *c, size_t idx);

qbid circ_sysqb(const struct circ *c, size_t idx);

qbid circ_ancqb(const struct circ *c, size_t idx);


#endif //PHASE2_CIRC_H
