# phase2 testing

Single-page reference for the test subsystem.  Sources:
`test/`.  Harness header: `test/test.h`.  Runner:
`test/run.c`.  Make integration: the `Testing`
section of `Makefile`.

## 1. Quick reference

Build and run the full suite:

```
make check
```

Subset:

```
make check-data            # all t-data_* tests
make check-paulis          # just t-paulis
make check-circ            # all t-circ*
```

Under MPI:

```
make check-mpi MPIRANKS=4  # wrap each test in mpirun -n 4
```

Slow tests (built separately, not part of default
`check`):

```
make check-slow
```

Sanitiser / valgrind variants:

```
make test-asan
make test-valgrind
make test-mpi-asan
```


## 2. Layout

Sources live under `test/`; compiled binaries land
under `build/test/` (mirrors the source tree, kept
out of git).

```
test/                  (sources)
  test.h               -- macro core (TEST_FAIL, TEST_ASSERT, ...)
  t-<subsys>[_<aspect>].c   -- one binary per logical unit
  t-data.h             -- per-subsystem fixture helpers (model)
  t-ref-coeff_matrix.py     -- Python cross-validator
  run.c                -- the parallel runner
  t-run.c              -- self-tests for the runner
  run-fixtures/*.c     -- synthetic fixtures driven by t-run
  data/                -- on-disk HDF5 fixtures for the C tests
  ref/                 -- in-tree Python reference oracles

build/test/            (binaries -- emitted by the build)
  t-<subsys>[_<aspect>]      -- test binary
  run                  -- parallel runner
  run-fixtures/{pass,fail,sleep,abort,banner}   -- t-run fixtures
```

C binaries follow the prefix convention
`t-<subsys>[_<aspect>].c`.  Examples:
`t-paulis.c`, `t-data_hamil.c`, `t-circ_trott_coeff.c`.
The matching binary lands at
`build/test/<same-name>` and is listed in the
Makefile's `TESTS` variable.

The slow test (`t-state_prep_coeff_large`) is declared
in `TESTS_SLOW` and is excluded from default `check`.

The Makefile target `check-tests-coverage` enforces
that every `test/t-*.c` source appears in TESTS or
TESTS_SLOW; it is wired as a prerequisite of
`build-test`, so a forgotten declaration fails the
build immediately.


## 3. `test.h` macros

```c
#include "test.h"

TEST_FAIL("unreachable branch hit (n=%zu)", n);
TEST_ASSERT(p != NULL, "alloc returned NULL for %zu bytes", sz);
TEST_EQ(rc, 0);
TEST_NEAR(observed, expected, 1e-12);
TEST_CNEAR(amp, 1.0 + 0.5*I, 1e-15);
```

| Macro                       | Use                                |
| --------------------------- | ---------------------------------- |
| `TEST_FAIL(fmt, ...)`       | unconditional abort                |
| `TEST_ASSERT(expr, fmt, ...)` | abort iff `expr` is false        |
| `TEST_EQ(a, b)`             | shorthand for `a == b`             |
| `TEST_NEAR(a, b, eps)`      | `fabs(a - b) <= eps` (double)      |
| `TEST_CNEAR(a, b, eps)`     | `cabs(a - b) <= eps` (_Complex)    |

All five are statements (`do { ... } while (0)`) -- do
not use in expression position.  Failure emits one
line to stderr in the GCC-error format:

```
<file>:<line>: FAIL [<func>]: <msg>
```

so editors recognise the prefix and jump directly to
the failing source location.  After the line, the
process exits with non-zero status.  No global cleanup
runs; the OS reaps the world.


## 4. Per-subsystem fixture headers

`test.h` itself stays macro-only.  Shared fixture
arrays, helper functions, and the
`test_fixture_copy`-style scaffolding go in a
per-subsystem header.  `test/t-data.h` is the
canonical model:

```c
/* test/t-data.h */
static const struct test_data {
        const char *path;
        size_t      n_orb, n_alpha, n_beta;
        ...
} TEST_DATA[] = { ... };

static void test_fixture_copy(const char *src, char *dst);
```

A test then `#include "t-data.h"` and iterates
`TEST_DATA[]`.  This keeps the macro core small and
keeps fixture wiring out of `test.h`.

The Makefile already lists `t-data.h` as a per-test
prerequisite, so a touch to the header rebuilds every
binary that consumes it.

Every test calls `world_init(NULL, NULL, SEED)` at
the top of `main` and `world_free()` before return.
Any `TEST_*` macro aborts via `exit`, leaving the
world uninitialised on the way out -- that is
intentional; the kernel reaps it.


## 5. The runner

The runner source is `test/run.c`; the build emits
the binary at `build/test/run`.  Standalone
single-file C program (libc + POSIX only, no
phase2 / MPI / HDF5 deps).  It takes a list of
test paths and runs them in parallel, capturing
per-test output, and prints a cargo-style summary.

```
./build/test/run [OPTIONS] -- TEST...
```

| Flag             | Default      | Effect                              |
| ---------------- | ------------ | ----------------------------------- |
| `--jobs=N`       | `nproc`      | parallel job count                  |
| `--mpiranks=N`   | 0 (off)      | wrap each test in `mpirun -n N`     |
| `--filter=GLOB`  | none         | fnmatch(3) against effective name   |
| `--verbose`, `-v`| off          | stream per-test stdout; forces `--jobs=1` |
| `--no-color`     | TTY-detected | force plain output                  |
| `--help`, `-h`   |              | usage + exit                        |

The *effective name* is the basename of the test path
with the leading `t-` and any `.py` suffix stripped:
`./build/test/t-data_mpi -> "data_mpi"`,
`./test/t-ref-coeff_matrix.py -> "ref-coeff_matrix"`.
Filter examples: `--filter='data*'`,
`--filter='*coeff*'`, `--filter='paulis'`.

### 5.1 Output

```
running 31 tests

test t-bitstring_index ............ ok (0.038s)
test t-circ_cache ................. ok (0.041s)
test t-data_mpi ................... FAILED (1.234s)
...

failures:

---- t-data_mpi stdout ----
<captured>
---- t-data_mpi stderr ----
test/t-data_mpi.c:108: FAIL [t_bcast_eq]: rank 1 diverges
(exit 1)

failures:
    t-data_mpi

test result: FAILED. 30 passed; 1 failed; finished in 6.32s
```

Per-test elapsed reported at millisecond precision so
sub-second variation between parallel tests is
visible.  Colour is on iff stdout is a TTY or
`FORCE_COLOR=1` is set; `--no-color` and pipe
redirection turn it off.  The plain format is
byte-identical apart from the ANSI escapes.

### 5.2 Exit codes

| Code | Meaning                                       |
| ---- | --------------------------------------------- |
| 0    | every dispatched test passed (or empty suite) |
| 1    | one or more tests failed                      |
| 2    | usage error (no `--` / no test paths)         |

Children killed by a signal are reported as
`128 + signum` in the failure dump; the runner itself
still exits 1.

### 5.3 Process model

Each test runs as a `fork(2)` + `execvp(3)` child.
stdout and stderr are captured to per-test tempfiles
under a `mkdtemp(3)` directory in `/tmp`; the captures
are surfaced under labelled headers in the failure
dump.  The parent uses `waitpid(-1, ...)` to reap
completions in arrival order and dispatches a new
child into each freed slot.  Children inherit the
parent's process group, so a `SIGINT` at the terminal
takes the whole run down cleanly without a custom
handler.

### 5.4 HDF5 file locking

HDF5 1.10+ takes an advisory `flock(2)` on every open;
parallel readers of committed fixtures intermittently
lose the race.  The runner sets
`HDF5_USE_FILE_LOCKING=FALSE` at startup, propagated
to every child via the environment.  All in-tree
tests open fixtures read-only, so suppressing the
lock is safe.

The architectural detail (process model, fd plumbing,
mkstemp choreography, Python wrapping, MPI wrapping)
is documented at the top of `test/run.c`.


## 6. Make integration

| Target                  | Effect                                  |
| ----------------------- | --------------------------------------- |
| `build-test`            | build all test binaries + `test/run`    |
| `build-test-slow`       | also build TESTS_SLOW                   |
| `check`                 | runner on TESTS + the .py harness       |
| `check-mpi`             | runner with `--mpiranks=$(MPIRANKS)`    |
| `check-slow`            | runner on TESTS_SLOW                    |
| `check-<filter>`        | runner with `--filter='<filter>*'`      |
| `check-tests-coverage`  | guard: every `t-*.c` is declared        |
| `test-asan`             | rebuild + run under ASan/UBSan          |
| `test-valgrind`         | rebuild + run under valgrind            |
| `test-mpi-asan`         | rebuild + run a small MPI subset under ASan |

The runner is invoked once per `make check`-family
target; it does its own parallel dispatch.  Make's
`-j` is therefore orthogonal -- not load-amplifying,
because the runner already saturates the cores.

`check-<filter>` is a pattern target; the stem is
threaded straight into the runner's `--filter` as
`<stem>*`.  Examples:

```
make check-data    -> --filter='data*'
make check-paulis  -> --filter='paulis*'
make check-circ    -> --filter='circ*'
```

Explicit targets (`check-mpi`, `check-slow`,
`check-tests-coverage`) win over this pattern.


## 7. Adding a new test

1. Create `test/t-<subsys>[_<aspect>].c` using the
   `test.h` macros.  If the test needs shared
   fixtures, include the per-subsystem fixture
   header (`test/t-data.h` or analogue).
2. Add the binary to `TESTS` in `Makefile`.  If it
   is a long-running test, add it to `TESTS_SLOW`
   instead.
3. `make check`.  The `check-tests-coverage`
   prerequisite of `build-test` will fail the build
   if step 2 was forgotten.

The Makefile already arranges:

- `-DDEBUG -g -Og` for every test binary,
- `-UDEBUG` for `t-log_release` (verifies the
  release-build strip of trace/debug macros),
- a generic dependency on `test/test.h`, `test/t-data.h`,
  the phase2/lib/ph2run object sets.


## 8. Adding a runner self-test

The runner has its own meta-test at `test/t-run.c`
which exercises `./test/run` via `popen(3)` against
synthetic fixtures under `test/run-fixtures/`
(pass, fail, sleep, abort, banner).  To add a new
scenario:

1. Drop a tiny C program under `test/run-fixtures/`
   that produces the controlled outcome you want
   (specific exit code, output pattern, runtime).
2. Add it to the `RUNFIX` Makefile variable so it
   gets built.
3. Add a scenario function in `test/t-run.c` that
   shells out to `./test/run` and asserts via
   `must_contain` / `must_not_contain` / `TEST_EQ`.
   Call it from `main`.

Self-tests are skipped automatically when the
runner is invoked under `mpirun` -- the meta-test
is single-process by design.


## 9. Sanitiser and valgrind targets

`test-asan` rebuilds the suite with
`-fsanitize=address,undefined -fno-omit-frame-pointer
-g` on both `EXTRA_CFLAGS` and `EXTRA_LDFLAGS`, then
runs `make check`.  OpenMPI's startup path triggers
internal "leaks" via libevent / libopen-pal that we
cannot fix from here; the rule disables leak
detection (`detect_leaks=0`) but keeps every other
ASan + UBSan diagnostic active.

`test-valgrind` runs the existing test binaries
under `valgrind --leak-check=full --track-origins=yes`
in a serial loop.  Does *not* go through the runner
(valgrind already serialises and produces its own
captures).

`test-mpi-asan` is a narrow build: rebuild with
ASan flags, then run a single representative MPI
test (`t-circ_trott_coeff`) under
`mpirun -n 4 -x ASAN_OPTIONS=...`.

All three targets are independent of the cargo-style
runner and reuse `$(TESTS)` directly.
