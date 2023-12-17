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

/**
 * Initialize global circuit environment.
 *
 * Although this function is exposed in the public API, the user doesn't need
 * to call it directly.  It will be called by circ_create() as well.
 *
 * Similarly, the function circ_shutdown() belows will be scheduled with
 * atexit() to clear the global environment during program exit.
 *
 * The initialization may take some time, especially when MPI mode is enabled.
 * It the user cannot wait longer than usual during the first call
 * to circ_init(), it is recomended that that circ_initialize() be called
 * directly at the beginning of the program.
 *
 * This function can be called multiple times from different threads.
 *
 * Return value:	 0	if the environment was successfully initialized
 *				by this function call
 *			 1	if the environment has already been initialized
 *				by another call to this function
 *			-1	in case of failure
 */
int
circ_initialize(void);

void
circ_shutdown(void);

struct circ *
circ_create(struct circuit *ct, void *data);

void
circ_destroy(struct circ *c);

void *
circ_data(const struct circ *c);

int
circ_report(struct circ const *c);

int
circ_reset(struct circ *c);

int
circ_run(struct circ *c);

size_t
circ_num_meaqb(const struct circ *c);

size_t
circ_num_sysqb(const struct circ *c);

size_t
circ_num_ancqb(const struct circ *c);

qbid
circ_meaqb(const struct circ *c, size_t idx);

qbid
circ_sysqb(const struct circ *c, size_t idx);

qbid
circ_ancqb(const struct circ *c, size_t idx);

#endif // PHASE2_CIRC_H
