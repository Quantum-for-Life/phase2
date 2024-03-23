#include <complex.h>

#include "qreg.h"

#include "circ.h"
#include "circ_private.h"

double circ_ops_prob0(struct circ *c, qbid qb)
{
	return 0.0; // calcProbOfOutcome(c->quest_qureg, qb, 0);
}

void circ_ops_blankstate(struct circ *c)
{
	qreg_zerostate(&c->quest_qureg);
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

void circ_ops_crotpauli(struct circ *c, int *paulis, double angle)
{
	struct paulis code_lo = paulis_new();
	for (u32 i = 0; i < c->quest_qureg.qb_lo; i++)
		paulis_set(&code_lo, paulis[i], i);
	struct paulis code_hi = paulis_new();
	for (u32 i = 0; i < c->quest_qureg.qb_hi; i++)
		paulis_set(&code_hi, paulis[i + c->quest_qureg.qb_lo], i);

	qreg_paulirot(&c->quest_qureg, code_hi, &code_lo, &angle, 1);
}
