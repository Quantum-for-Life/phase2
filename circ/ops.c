#include <complex.h>

#include "types.h"

#include "qreg.h"

#include "circ.h"
#include "circ_private.h"

double circ_ops_prob0(struct circ *c, qbid qb)
{
	return 0.0; // calcProbOfOutcome(c->quest_qureg, qb, 0);
}

void circ_ops_blankstate(struct circ *c)
{
	qreg_blankstate(&c->quest_qureg);
}

void circ_ops_setsysamp(struct circ *c, size_t idx, _Complex double amp)
{
	double	  amps[2]   = { creal(amp), cimag(amp) };
	long long start_ind = idx << circ_num_meaqb(c);
	qreg_setamp(&c->quest_qureg, start_ind, amps);
}

void circ_ops_getsysamp(struct circ *c, size_t idx, _Complex double *amp)
{
	double	  amps[2];
	long long start_ind = idx << circ_num_meaqb(c);
	qreg_getamp(&c->quest_qureg, start_ind, &amps);

	*amp = amps[0] + _Complex_I * amps[1];
}

void circ_ops_multirotpauli(struct circ *c, struct paulis code_hi,
	const struct paulis *codes_lo, const fl *angles, size_t num_codes)
{
	qreg_paulirot(&c->quest_qureg, code_hi, codes_lo, angles, num_codes);
}

void circ_ops_rotpauli(struct circ *c, struct paulis code, fl angle)
{
	struct paulis code_hi, code_lo;
	paulis_split(code, c->quest_qureg.qb_lo, c->quest_qureg.qb_hi, &code_lo,
		&code_hi);
	paulis_shr(&code_hi, c->quest_qureg.qb_lo);

	circ_ops_multirotpauli(c, code_hi, &code_lo, &angle, 1);
}
