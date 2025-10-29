#ifndef CIRC_H
#define CIRC_H

#include <stddef.h>
#include <stdint.h>

#include "phase2/data.h"
#include "phase2/paulis.h"
#include "phase2/qreg.h"

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

struct circ_values {
	_Complex double *z;
	size_t len;
};

struct circ {
	struct circ_hamil hm;
	struct circ_muldet md;
	struct circ_cache *cache;
	struct circ_values vals;
	struct qreg reg;
};

int circ_hamil_init(struct circ_hamil *hm, uint32_t qb, size_t len);
void circ_hamil_free(struct circ_hamil *hm);
void circ_hamil_sort_lex(struct circ_hamil *hm);

int circ_muldet_init(struct circ_muldet *md, size_t len);
void circ_muldet_free(struct circ_muldet *md);

void circ_prog_init(struct circ_prog *prog, size_t len);
void circ_prog_tick(struct circ_prog *prog);

int circ_values_init(struct circ_values *vals, size_t len);
void circ_values_free(struct circ_values *vals);

int circ_init(struct circ *ct, data_id fid, size_t vals_len);
void circ_free(struct circ *ct);
int circ_prepst(struct circ *ct);
int circ_step(struct circ *ct, const struct circ_hamil *hm, double omega);
int circ_step_reverse(
	struct circ *ct, const struct circ_hamil *hm, double omega);
_Complex double circ_measure(struct circ *ct);

#endif // CIRC_H
