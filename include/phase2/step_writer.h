#ifndef PHASE2_STEP_WRITER_H
#define PHASE2_STEP_WRITER_H

#include <stddef.h>
#include <complex.h>

/*
 * Per-step result-writer callback used by the
 * trott / trott2 / qdrift / cmpsit algorithms.
 *
 * `write` is invoked once per simulation step on rank 0
 * (followers receive a no-op or short-circuit through the
 * implementation).  `ctx` is opaque to phase2; the caller
 * (typically ph2run) supplies it together with the
 * function pointer.  Returns 0 on success, -1 on error.
 *
 * Callers may pass a NULL `struct phase2_step_writer *`
 * to an algorithm `*_init` to disable per-step output;
 * the algorithm then runs entirely in memory.
 */
struct phase2_step_writer {
	void *ctx;
	int (*write)(void *ctx, size_t step_idx, _Complex double z);
};

#endif /* PHASE2_STEP_WRITER_H */
