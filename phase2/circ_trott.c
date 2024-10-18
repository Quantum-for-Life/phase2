#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "phase2/world.h"

#include "circ_impl.h"

#define MAX_CACHE_CODES (1024)


struct circ_trott {
	size_t nqb;
	struct qreg reg;

	const struct circ_trott_data *data;

	double prod[2];

	struct code_cache {
		struct paulis code_hi;
		struct paulis codes_lo[MAX_CACHE_CODES];
		double angles[MAX_CACHE_CODES];
		size_t ncodes;
	} cache;
};

static int circ_create(struct circ_trott *c, const struct circ_trott_data *data,
	const size_t nqb)
{
	struct qreg reg;
	if (qreg_init(&reg, nqb) < 0)
		return -1;

	c->nqb = nqb;
	c->data = data;
	c->reg = reg;

	return 0;
}

static void circ_destroy(struct circ_trott *c)
{
	qreg_destroy(&c->reg);
}

int circ_trott_data_init(struct circ_trott_data *cd, size_t ntsteps)
{
	cd->ntsteps = ntsteps;
	cd->tsteps[0] = malloc(sizeof(double) * 2 * ntsteps);
	if (cd->tsteps[0] == NULL)
		return -1;
	cd->tsteps[1] = cd->tsteps[0] + ntsteps;

	return 0;
}

void circ_trott_data_destroy(struct circ_trott_data *cd)
{
	circ_multidet_destroy(&cd->multidet);
	circ_hamil_destroy(&cd->hamil);

	free(cd->tsteps[0]);
}

int circ_trott_data_init_from_file(
	struct circ_trott_data *cd, size_t ntsteps, data_id fid)
{
	if (circ_trott_data_init(cd, ntsteps) < 0)
		return -1;

	int rc = circ_hamil_init_from_file(&cd->hamil, fid);
	rc |= circ_multidet_init_from_file(&cd->multidet, fid);
	data_circ_trott_getttrs(fid, &cd->delta);

	return rc;
}

static int circuit_prepst(struct circ_trott *c)
{
	const struct circ_multidet *md = &c->data->multidet;

	qreg_zero(&c->reg);
	for (size_t i = 0; i < md->ndets; i++) {
		const _Complex double cf =
			md->dets[i].cf[0] + _Complex_I * md->dets[i].cf[1];
		qreg_setamp(&c->reg, md->dets[i].idx, cf);
	}

	return 0;
}

static void trott_step(struct circ_trott *c, const double omega)
{
	const struct circ_hamil *hamil = &c->data->hamil;

	struct code_cache cache = c->cache;
	cache.ncodes = 0;

	for (size_t i = 0; i < hamil->nterms; i++) {
		const double angle = omega * hamil->cfs[i];
		const struct paulis code = hamil->ops[i];

		struct paulis code_hi, code_lo;
		paulis_split(
			code, c->reg.nqb_lo, c->reg.nqb_hi, &code_lo, &code_hi);

		if (cache.ncodes == 0) {
			cache.code_hi = code_hi;
			cache.codes_lo[0] = code_lo;
			cache.angles[0] = angle;
			cache.ncodes++;
			continue;
		}

		if (paulis_eq(cache.code_hi, code_hi) &&
			cache.ncodes < MAX_CACHE_CODES) {
			const size_t k = cache.ncodes++;
			cache.codes_lo[k] = code_lo;
			cache.angles[k] = angle;
			continue;
		}

		log_trace(
			"paulirot, term: %zu, num_codes: %zu", i, cache.ncodes);
		qreg_paulirot(&c->reg, cache.code_hi, cache.codes_lo,
			cache.angles, cache.ncodes);

		cache.ncodes = 1;
		cache.code_hi = code_hi;
		cache.codes_lo[0] = code_lo;
		cache.angles[0] = angle;
	}

	log_trace("paulirot, last term group, num_codes: %zu", cache.ncodes);

	if (cache.ncodes > 0)
		qreg_paulirot(&c->reg, cache.code_hi, cache.codes_lo,
			cache.angles, cache.ncodes);
}

static int circ_effect(struct circ_trott *c)
{
	const double t = c->data->delta;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;

	trott_step(c, t);

	return 0;
}

static int circ_measure(struct circ_trott *c)
{
	const struct circ_multidet *md = &c->data->multidet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->ndets; i++) {
		_Complex double a;
		qreg_getamp(&c->reg, md->dets[i].idx, &a);

		const _Complex double damp =
			md->dets[i].cf[0] + _Complex_I * md->dets[i].cf[1];
		pr += a * conj(damp);
	}
	c->prod[0] = creal(pr);
	c->prod[1] = cimag(pr);

	return 0;
}

int circ_simulate(struct circ *c)
{
	int rt = -1;

	size_t prog_percent = 0;
	const size_t nqb = cd->hamil.nqb;

	struct circ_trott c;
	if (circ_create(&c, cd, nqb) < 0)
		goto exit_circ_create;
	circuit_prepst(&c);

	for (size_t i = 0; i < cd->ntsteps; i++) {
		size_t percent = i * 100 / cd->ntsteps;
		if (percent > prog_percent) {
			prog_percent = percent;
			log_info("Progress: %zu\% (trott_step: %zu)", percent,
				i);
		}

		if (circ_effect(&c) < 0)
			goto exit_circ_effect;
		circ_measure(&c);
		cd->tsteps[0][i] = c.prod[0];
		cd->tsteps[1][i] = c.prod[1];
	}

	rt = 0; /* Success. */

exit_circ_effect:
	circ_destroy(&c);
exit_circ_create:

	return rt;
}
