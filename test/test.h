#ifndef TEST_H
#define TEST_H

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Minimal C test harness.
 *
 * Three core macros:
 *   TEST_FAIL(fmt, ...)        unconditional abort with
 *                              file:line:func and message;
 *                              exits with non-zero status.
 *   TEST_ASSERT(expr, fmt, ...) abort with the formatted
 *                              message iff `expr` is false.
 *   TEST_EQ(a, b)              shorthand for the common
 *                              `expected a == b` check.
 *
 * All three macros are statements -- never use in
 * expression position.  __VA_OPT__ keeps the comma off when
 * the variadic list is empty (C23 / GNU C).
 *
 * On failure each macro emits one line to stderr in the
 * GCC-error format:
 *
 *   <file>:<line>: FAIL [<func>]: <msg>
 *
 * so editors recognise the leading `file:line:` and jump
 * directly to the failing source location.
 *
 * Conventions for new tests:
 *
 *   - Binaries live in `test/` and are named
 *     `t-<subsys>[_<aspect>].c` (e.g. `t-data_hamil.c`,
 *     `t-circ_trott_coeff.c`); the matching binary lands
 *     under the same name, listed in the Makefile's TESTS
 *     variable.
 *   - Shared fixture helpers go in a per-subsystem header
 *     (`test/t-data.h` is the canonical model -- it carries
 *     the `TEST_DATA[]` fixture array and the
 *     `test_fixture_copy` helper).  `test.h` itself stays
 *     macro-only.
 *   - Each test calls `world_init(nullptr, nullptr, SEED)`
 *     at the top of `main` and `world_free()` before
 *     return.  Failure path: any TEST_* macro aborts via
 *     `exit`, so the world is left uninitialised on the
 *     way out -- that is intentional, the OS reaps it.
 */

#define TEST_FAIL(fmt, ...)                                                    \
	do {                                                                   \
		fprintf(stderr, "%s:%d: FAIL [%s]: " fmt "\n", __FILE__,       \
			__LINE__, __func__ __VA_OPT__(, ) __VA_ARGS__);        \
		exit(-1);                                                      \
	} while (0)

#define TEST_ASSERT(expr, fmt, ...)                                            \
	do {                                                                   \
		if (!(expr))                                                   \
			TEST_FAIL(fmt __VA_OPT__(, ) __VA_ARGS__);             \
	} while (0)

#define TEST_STR(a) #a
#define TEST_XSTR(a) TEST_STR(a)

#define TEST_EQ(a, b)                                                          \
	do {                                                                   \
		TEST_ASSERT((a) == (b), "expected %s == %s", TEST_XSTR(a),     \
			TEST_XSTR(b));                                         \
	} while (0)

/*
 * Tolerance-based equality for doubles.  Aborts when
 * |a - b| > eps; the failure message carries the
 * stringified expressions plus their observed values
 * and the violated diff.
 */
#define TEST_NEAR(a, b, eps)                                                   \
	do {                                                                   \
		const double _a = (a), _b = (b), _eps = (eps);                 \
		if (fabs(_a - _b) > _eps)                                      \
			TEST_FAIL("|%s - %s| = %g > %g (%g vs %g)",            \
				TEST_XSTR(a), TEST_XSTR(b), fabs(_a - _b),     \
				_eps, _a, _b);                                 \
	} while (0)

/*
 * Tolerance-based equality for _Complex doubles.  Same
 * shape as TEST_NEAR but uses cabs and prints both
 * components on failure.
 */
#define TEST_CNEAR(a, b, eps)                                                  \
	do {                                                                   \
		const _Complex double _a = (a), _b = (b);                      \
		const double _eps = (eps);                                     \
		if (cabs(_a - _b) > _eps)                                      \
			TEST_FAIL("|%s - %s| = %g > %g"                        \
				  " ((%g+%gi) vs (%g+%gi))",                   \
				TEST_XSTR(a), TEST_XSTR(b), cabs(_a - _b),     \
				_eps, creal(_a), cimag(_a), creal(_b),         \
				cimag(_b));                                    \
	} while (0)

#endif /* TEST_H */
