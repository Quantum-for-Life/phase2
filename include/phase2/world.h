#ifndef WORLD_H
#define WORLD_H

#include <stdint.h>

#include "log.h"
#include "xoshiro256ss.h"

enum world_stat {
	WORLD_UNDEF = -1,
	WORLD_READY = 0,
	WORLD_DONE = 1,
	WORLD_ERR = 2,
};

#ifndef PHASE2_BACKEND
#define PHASE2_BAKCEND (0)
#endif

#if PHASE2_BACKEND == 0
#define WORLD_BACKEND "qreg"
#elif PHASE2_BACKEND == 1
#define WORLD_BACKEND "QuEST"
#elif PHASE2_BACKEND == 2
#define WORLD_BACKEND "CUDA"
#endif /* PHASE2_BACKEND */

/*
 * Global world info stored as static data inside the module.
 * Users can keep a local copy of this when populated by world_info().
 */
struct world {
	enum world_stat stat;
	int size;
	int rank;

	uint64_t seed;
	struct xoshiro256ss rng;

	void *data; /* opaque handle to alternative engine environments */
};

/* Initialize the world with command line parameters and a seed for PRNG.
 *
 * This should be called exactly once at the begging of the program.  The MPI
 * world communicator is initialized here as well as the logging facility.
 *
 * The PRNG is initialized with the seed.  The seed must not be zero.  The state
 * of the PRNG is deterministically split among MPI processes, so that each
 * process has it own state.
 *
 * Returns:
 * 	WORLD_READY	- in case of success
 *	WORLD_ERR	- if an error occurred, e.g. seed is zero
 */
int world_init(int *argc, char ***argv, uint64_t seed);

/* Destroy the world (global simulation environment).
 *
 * This should be called exactly once at the end of the program.
 *
 * This function deinitialized the log facility as well, so no log messages
 * will be recorded after calling this function.
 *
 * Returns:
 *	WORLD_DONE	- in case of success
 *	WORLD_ERR	- in case of error
 */
int world_free(void);

/* Get information about the world.
 *
 * This populates the supplied struct with the information about the global
 * static world structure.  It does not change the world state or parameters.
 *
 * Returns:
 *	The same value as stored in wd->stat after the call.
 */
int world_info(struct world *wd);

int world_backend_init(struct world *wd);

void world_backend_destroy(struct world *wd);

#endif /* WORLD_H */
