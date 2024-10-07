#ifndef WORLD_H
#define WORLD_H

#include <stdint.h>

#include "xoshiro256ss.h"

#define PHASE2_LOG_ENVVAR "PHASE2_LOG"

enum world_log_level {
	LOG_TRACE,
	LOG_DEBUG,
	LOG_INFO,
	LOG_WARN,
	LOG_ERROR,
	LOG_FATAL
};

/* Call world_init() before using this. */
void world_log(int level, const char *fmt, ...);

#define log_trace(...) world_log(LOG_TRACE, __VA_ARGS__)
#define log_debug(...) world_log(LOG_DEBUG, __VA_ARGS__)
#define log_info(...) world_log(LOG_INFO, __VA_ARGS__)
#define log_warn(...) world_log(LOG_WARN, __VA_ARGS__)
#define log_error(...) world_log(LOG_ERROR, __VA_ARGS__)
#define log_fatal(...) world_log(LOG_FATAL, __VA_ARGS__)

enum world_stat {
	WORLD_UNDEF	= -1,
	WORLD_READY	=  0,
	WORLD_DONE	=  1,
	WORLD_ERR	=  2,
};

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

	void *data;	/* opaque handle to alternative engine environments */
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
int world_destroy(void);

/* Get information about the world.
 *
 * This populates the supplied struct with the information about the global
 * static world structure.  It does not change the world state or parameters.
 *
 * Returns:
 *	The same value as stored in wd->stat after the call.
 */
int world_info(struct world *wd);

#endif /* WORLD_H */
