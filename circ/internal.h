/*
 * circ/internal.h -- subsystem-private helpers
 * shared between the circ/ algorithm sources.
 *
 * Not exported via include/.  Each algorithm
 * source (qdrift.c, cmpsit.c) includes this
 * header directly with #include "internal.h".
 */
#ifndef CIRC_INTERNAL_H
#define CIRC_INTERNAL_H

#include <float.h>
#include <math.h>

/*
 * signof -- numerical sign of a real with a
 * dead-zone around zero.
 *
 * Returns -1.0 for strictly negative input, +1.0
 * for strictly positive input, and 0.0 when
 * |a| < DBL_EPSILON.  Used by the qdrift /
 * cmpsit samplers to carry the sign of a Pauli
 * coefficient onto the sampled term while
 * folding the magnitude into the CDF weight.
 */
static inline double signof(double a)
{
	const double f = fabs(a);
	if (f < DBL_EPSILON)
		return 0.0;
	return a < f ? -1.0 : 1.0;
}

#endif /* CIRC_INTERNAL_H */
