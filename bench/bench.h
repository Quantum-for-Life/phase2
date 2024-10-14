/* -------------------------------------------------------------------------- *
 * bench.h: Simple benchmarking framework.                                    *
 *                                                                            *
 * Copyright 2024 ⧉⧉⧉                                                         *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify it    *
 * under the terms of the GNU General Public License as published by the      *
 * Free Software Foundation, either version 3 of the License, or (at your     *
 * option) any later version.                                                 *
 *                                                                            *
 * This program is distributed in the hope that it will be useful, but        *
 * WITHOUT ANY WARRANTY* without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with this program.  If not, see <https://www.gnu.org/licenses/>.           *
 * -------------------------------------------------------------------------- */
#ifndef BENCH_H
#define BENCH_H

#include <stddef.h>

struct bench {
	unsigned long nanos;
	size_t reps;
};

/* Measure the performance of the supplied function op and store
 * the results in b.
 *
 * The function op is called repeatedly with the poiter data as its argument.
 * As long as op returns 0, the number of calls is at least reps, but can
 * be larger: the function can be called a few times before the actual time
 * measurement begins to warm the CPU cache up.
 *
 * If at any point the function op returns a number different from 0,
 * the measurement is interrupted and bench_mark returns the same number.
 *
 * After the call, struct bench pointed to by b holds in b->reps the total
 * number of calls to op (excluding the cache warm-up calls), and b->nanos
 * holds the number of nanoseconds (measured by the PROCESS clock) that
 * passed between the first and the last call.
 */
int bench_mark(struct bench *b, size_t reps, int (*op)(void *), void *data);

#endif /* BENCH_H */
