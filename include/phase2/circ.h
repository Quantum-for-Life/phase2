#ifndef CIRC_H
#define CIRC_H

#include <stddef.h>
#include <stdint.h>

#include "phase2/data.h"
#include "phase2/paulis.h"
#include "phase2/qreg.h"

#define CIRC_CACHE_CODES_MAX UINT64_C(0x0400)

struct circ_cache {
	uint32_t qb_lo, qb_hi;
	struct paulis *codes_lo, code_hi;
	double *phis;
	size_t len;
};

struct circ_hamil {
	uint32_t qb;
	struct circ_hamil_term {
		double cf;
		struct paulis op;
	} *terms;
	size_t len;
};

struct circ_muldet {
	struct {
		uint64_t idx;
		_Complex double cf;
	} *dets;
	size_t len;
};

struct circ_prog {
	unsigned pc; /* percent */
	size_t i, len;
};

struct circ {
	struct circ_hamil hm;
	struct circ_muldet md;
	struct circ_cache cache;
	struct qreg reg;
};

int circ_cache_init(struct circ_cache *ch, uint32_t qb_lo, uint32_t qb_hi);
void circ_cache_free(struct circ_cache *ch);
int circ_cache_insert(struct circ_cache *ch, struct paulis code, double phi);
void circ_cache_flush(struct circ_cache *ch,
	void (*op)(struct paulis, struct paulis *, double *, size_t, void *),
	void *data);

int circ_hamil_init(struct circ_hamil *hm, uint32_t qb, size_t len);
void circ_hamil_free(struct circ_hamil *hm);
void circ_hamil_sort_lex(struct circ_hamil *hm);

int circ_muldet_init(struct circ_muldet *md, size_t len);
void circ_muldet_free(struct circ_muldet *md);

void circ_prog_init(struct circ_prog *prog, size_t len);
void circ_prog_tick(struct circ_prog *prog);

int circ_init(struct circ *ct, data_id fid);
void circ_free(struct circ *ct);
int circ_prepst(struct circ *ct);
int circ_step(struct circ *ct, const struct circ_hamil *hm, double omega);
_Complex double circ_measure(struct circ *ct);

#endif // CIRC_H
