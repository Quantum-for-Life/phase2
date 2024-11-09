#ifndef CIRC_H
#define CIRC_H

#include <stddef.h>
#include <stdint.h>

#include "phase2/data.h"
#include "phase2/paulis.h"
#include "phase2/qreg.h"

#define MAX_CACHE_CODES UINT64_C(0x0400)

struct circ_cache {
	uint32_t qb_lo, qb_hi;
	struct paulis *codes_lo, code_hi;
	double *phis;
	size_t n;
};

int circ_cache_init(struct circ_cache *ch, size_t qb_lo, size_t qb_hi);

void circ_cache_destroy(struct circ_cache *ch);

int circ_cache_insert(struct circ_cache *ch, struct paulis code, double phi);

void circ_cache_flush(struct circ_cache *ch,
	void (*op)(struct paulis, struct paulis *, double *, size_t, void *),
	void *data);

struct circ_hamil {
	size_t nqb;

	struct circ_hamil_term {
		struct paulis op;
		double cf;
	} *terms;

	size_t nterms;
};

void circ_hamil_sort_lex(struct circ_hamil *hm);

struct circ_muldet {
	struct {
		uint64_t idx;
		_Complex double cf;
	} *dets;

	size_t ndets;
};

struct circ {
	struct circ_hamil hamil;
	struct circ_muldet muldet;
	struct circ_cache cache;
	struct qreg reg;

	int (*simul)(struct circ *);
};

int circ_init(struct circ *ct, data_id fid, int (*simul)(struct circ *));

void circ_destroy(struct circ *ct);

int circ_prepst(struct circ *ct);

int circ_step(struct circ *ct, const struct circ_hamil *hm, double omega);

_Complex double circ_measure(struct circ *ct);

int circ_simul(struct circ *ct);


#endif // CIRC_H
