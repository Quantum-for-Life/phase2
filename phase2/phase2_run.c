#include "c23_compat.h"
#include <complex.h>
#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"

#include "phase2/phase2_run.h"
#include "circ_cache.h"

#define PH2RUN_SEED UINT64_C(0xa7c3e1b9d4f20586)
#define PH2RUN_MAX_QB (30)

static int world_ready;

static void cleanup(void)
{
	if (world_ready)
		world_free();
}

static int ensure_world(void)
{
	struct world_info wd;

	world_info(&wd);
	if (wd.stat == WORLD_READY)
		return 0;

	if (world_init(nullptr, nullptr, PH2RUN_SEED) != WORLD_READY)
		return -1;

	world_ready = 1;
	atexit(cleanup);

	return 0;
}

static int parse_basis(const char *s, uint32_t nqb, uint64_t *out)
{
	*out = 0;
	for (uint32_t k = 0; k < nqb; k++) {
		switch (s[k]) {
		case '0':
			break;
		case '1':
			*out |= UINT64_C(1) << k;
			break;
		default:
			return -1;
		}
	}

	return 0;
}

static int parse_pauli(const char *s, uint32_t nqb, struct paulis *out)
{
	char buf[4096];
	size_t len = strlen(s);

	if (len >= sizeof(buf))
		return -1;

	memcpy(buf, s, len + 1);
	*out = paulis_new();

	char *tok = strtok(buf, " \t\n");
	while (tok) {
		int op;
		switch (tok[0]) {
		case 'I':
			op = PAULI_I;
			break;
		case 'X':
			op = PAULI_X;
			break;
		case 'Y':
			op = PAULI_Y;
			break;
		case 'Z':
			op = PAULI_Z;
			break;
		default:
			return -1;
		}

		char *end;
		unsigned long idx = strtoul(tok + 1, &end, 10);
		if (end == tok + 1 || *end != '\0')
			return -1;
		if (idx >= nqb)
			return -1;

		paulis_set(out, op, (uint32_t)idx);
		tok = strtok(nullptr, " \t\n");
	}

	return 0;
}

static void flush_cb(struct paulis code_hi,
	const struct paulis *codes_lo, double *phis,
	size_t ncodes, void *data)
{
	struct qreg *reg = data;

	qreg_paulirot(reg, code_hi, codes_lo, phis, ncodes);
}

int phase2_run(const char **pauli_strs,
	const double *coeffs, size_t nterms, double delta,
	const char *psi_str, double *out_re, double *out_im)
{
	uint32_t nqb;
	uint64_t idx;
	struct qreg reg;

	if (!psi_str)
		return -1;
	nqb = (uint32_t)strlen(psi_str);
	if (nqb == 0 || nqb > PH2RUN_MAX_QB)
		return -1;

	if (ensure_world() < 0)
		return -1;
	if (parse_basis(psi_str, nqb, &idx) < 0)
		return -1;
	if (qreg_init(&reg, nqb) < 0)
		return -1;

	qreg_zero(&reg);
	qreg_setamp(&reg, idx, 1.0);

	struct circ_cache *cache = circ_cache_init(reg.qb_hi, reg.qb_lo);
	if (!cache)
		goto err_cache;

	for (size_t k = 0; k < nterms; k++) {
		struct paulis code;

		if (parse_pauli(pauli_strs[k], nqb, &code) < 0)
			goto err_pauli;

		double phi = delta * coeffs[k];
		if (circ_cache_insert(cache, code, phi) != 0) {
			circ_cache_flush(cache, flush_cb, &reg);
			if (circ_cache_insert(cache, code, phi) < 0)
				goto err_pauli;
		}
	}
	circ_cache_flush(cache, flush_cb, &reg);

	_Complex double z;
	qreg_getamp(&reg, idx, &z);
	*out_re = creal(z);
	*out_im = cimag(z);

	circ_cache_free(cache);
	qreg_free(&reg);

	return 0;

err_pauli:
	circ_cache_free(cache);
err_cache:
	qreg_free(&reg);

	return -1;
}
