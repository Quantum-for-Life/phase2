#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "phase2/paulis.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "test.h"

static struct world WD;

#define WIDTH (64)
#define MARGIN (1.0e-10)

#define NUM_QUBITS (12)
#define NUM_AMPS (1UL << NUM_QUBITS)
static _Complex double AMPS[NUM_AMPS];

#define SEED UINT64_C(0xd3b9268b8737ddc0)
static struct xoshiro256ss RNG;

static double rand_double(void)
{
	uint64_t x = xoshiro256ss_next(&RNG);

	return (x >> 11) * 0x1.0p-53;
}

static enum pauli_op rand_pauli_op(void)
{
	enum pauli_op x = (int)(xoshiro256ss_next(&RNG) % 4);

	return x;
}

void TEST_MAIN(void)
{
	world_init((void *)0, (void *) 0);
	world_info(&WD);
	log_info("MPI world size: %d", WD.size);

	xoshiro256ss_init(&RNG, SEED);

	world_fin();
}
