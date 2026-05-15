# phase2 logging

Single-page reference for the logging subsystem.  Header:
`include/log.h`.  Implementation: `lib/log.c`.  Companion
tests: `test/t-log.c`, `test/t-log_release.c`.

## 1. Quick reference

```c
#define LOG_SUBSYS "myunit"
#include "log.h"

log_trace("hot-loop event id=%llu", id);  /* compiled in only with -DDEBUG */
log_debug("init: alpha=%u beta=%u", a, b);/* compiled in only with -DDEBUG */
log_info ("running simulation: %zu steps", n);
log_warn ("retrying after transient error %d", rc);
log_error("open(%s) failed: %s", path, strerror(errno));
log_fatal("unrecoverable: %s", reason);
```

Default output:

```
2026-05-13 14:22:11.034 INFO  [myunit] running simulation: 100 steps
```

Levels at or below INFO (`trace`/`debug`/`info`) go to
**stdout**; `warn`/`error`/`fatal` go to **stderr**.


## 2. Levels and policy

The level determines both verbosity and routing.  Pick the
right one from the table:

| Level  | Owner          | Log this                  | Do NOT log             |
| ------ | -------------- | ------------------------- | ---------------------- |
| trace  | phase2/, circ/ | batch / kernel entries    | per-amplitude state    |
| debug  | phase2/, circ/ | init, attribute reads     | already-handled events |
| info   | ph2run/, circ/ | options, paths, progress  | per-step engine chatter|
| warn   | anywhere       | recoverable degradation   | normal events          |
| error  | anywhere       | failing return path       | events already caught  |
| fatal  | anywhere       | about to `exit` / `abort` | anything recoverable   |

Full guidance on what counts as "batch boundary" vs
"per-amplitude" lives in §7 below.

The phase2 library reports mostly **below** info; ph2run
reports mostly **at** info.  Algorithm files in `circ/`
emit info-level progress and lifecycle banners.


## 3. Output format

```
[r=R] YYYY-MM-DD HH:MM:SS.mmm LEVEL [subsys] message
```

Field order, left to right:

- **rank prefix** `[r=R]`: only when `PHASE2_LOG_ALL=1`
  AND the calling rank is non-zero.  Rank 0 stays
  unprefixed in every mode, so it remains the dominant
  voice and is easy to scan with `grep -v '\[r='`.
- **timestamp**: local time, millisecond precision via
  `clock_gettime(CLOCK_REALTIME, ...)` +
  `localtime_r`.
- **level**: five-character right-padded label
  (`TRACE`/`DEBUG`/`INFO `/`WARN `/`ERROR`/`FATAL`).
- **subsystem tag**: per-file label, set by `#define
  LOG_SUBSYS "<name>"` before `#include "log.h"`.
  Defaults to `?` if missing; reviewers flag any `[?]`.
- **message**: printf-style; format is validated by GCC
  via `__attribute__((format(printf, …)))` on
  `log_emit`.


## 4. Environment variables

| Variable          | Default | Accepted values                          |
| ----------------- | ------- | ---------------------------------------- |
| `PHASE2_LOG`      | info    | `trace`/`debug`/`info`/`warn`/`error`/`fatal` |
| `PHASE2_LOG_ALL`  | unset   | any non-empty value enables rank > 0      |

`PHASE2_LOG` is case-insensitive and only the first
character is consulted, so `t`, `Trace`, and `TRACE`
parse the same.

Unparseable `PHASE2_LOG` values silently fall back to the
default — same behaviour as the prior release.


## 5. MPI semantics

Rank 0 emits unconditionally.  Other ranks emit only when
`PHASE2_LOG_ALL=1`.  When all-ranks mode is on, non-rank-0
lines begin with `[r=<rank>] ` so the operator can
`sort` the captured log to interleave correctly:

```
[r=0] 2026-05-13 14:22:11.034 INFO  [world] backend: qreg
[r=2] 2026-05-13 14:22:11.041 DEBUG [qreg]  paulirot batch n=64
[r=1] 2026-05-13 14:22:11.041 ERROR [qreg]  send timed out
```

The rank is read from `world_info()` on each emit, so the
flag and rank are always current.

**SLURM `.out` flushing.** `log_init()` calls
`setvbuf(stdout, NULL, _IOLBF, 0)` and
`setvbuf(stderr, NULL, _IONBF, 0)`, and the emit path
`fflush()`es after each line.  An operator running
`tail -f slurm-NNN.out` sees every log line as it is
produced, not in 4 KB chunks buffered until end of run.


## 6. Build-time gating

**Release builds strip `log_trace` and `log_debug` to
`((void)0)`.**  The header gates them on `DEBUG`:

```c
#ifdef DEBUG
#define log_trace(...) log_at(LOG_TRACE, __VA_ARGS__)
#define log_debug(...) log_at(LOG_DEBUG, __VA_ARGS__)
#else
#define log_trace(...) ((void)0)
#define log_debug(...) ((void)0)
#endif
```

Production binaries (`make build` / `make all`) do not
define `DEBUG`.  Setting `PHASE2_LOG=trace` on such a
binary has **no effect** on trace/debug emission — the
call sites have been compiled out, no load and no branch
remain.

To trace a production-class run, rebuild with debug
flags:

```sh
make distclean
make debug
```

Tests always build with `-DDEBUG -g -Og` (see Makefile);
the dedicated `test/t-log_release` binary cancels this
with a trailing `-UDEBUG` to verify the strip.


## 7. Policy for new code

These rules apply to every patch that adds or modifies a
function:

1. **Log non-trivial work.**  Every new function that does
   I/O, allocates a non-trivial resource, branches on
   configuration, or can fail **must** emit at least one
   log line.  Pick the level from §2.
2. **No silent error returns.**  Every error-return path
   emits `log_error` with the failing operation and the
   context an operator needs (path, index, return code).
   `return -1;` without a prior log_error is a review
   blocker.
3. **Subsystem at module entry.**  Every public entry
   point of a new module starts with one `log_debug`
   describing its arguments (engine layer) or
   `log_info` describing intent (driver layer).
4. **Tag your file.**  Every translation unit that calls
   `log_*` must `#define LOG_SUBSYS "<name>"` immediately
   before `#include "log.h"`.  A `[?]` in test capture is
   a review blocker.
5. **Never log per-amplitude.**  The release strip
   protects performance; this rule protects readability
   of debug captures.  Trace at the *outer* batch /
   kernel call, never inside the amplitude `for` loop.
6. **Format specifiers must match types.**  GCC enforces
   this via the `format(printf, …)` attribute on
   `log_emit`; use `<inttypes.h>` macros (`PRIu64`, etc.)
   for fixed-width integers.

This policy is enforced by review until a tooling check
(clang-tidy, pre-commit grep) lands.


## 8. Adding logging to existing code

Quick recipe when editing a file that does not yet emit:

1. Add `#define LOG_SUBSYS "<short-name>"` immediately
   before `#include "log.h"` (or `#include "phase2.h"`
   if that is what the file uses).  Pick the name from
   the existing convention: directory-aligned, lowercase,
   no spaces.
2. Identify the points that match §7 rule 1 (I/O,
   allocation, branches, failure).  Add the right level.
3. Recompile.  GCC will reject any format-spec mismatch.
4. Eyeball one run with `PHASE2_LOG=debug` to confirm the
   new lines have the right subsystem tag and level.


## 9. Examples

**Small run (info, default).**

```sh
./build/ph2run/ph2run -S test/data/H2O_CAS56.h5 trott -D 0.1 -s 4
```

```
2026-05-13 14:22:11.034 INFO  [world]  *** Init ***
2026-05-13 14:22:11.034 INFO  [world]  World size: 1
2026-05-13 14:22:11.034 INFO  [world]  Backend: qreg
2026-05-13 14:22:11.034 INFO  [ph2run] subcommand: trott
2026-05-13 14:22:11.034 INFO  [ph2run] simul file: test/data/H2O_CAS56.h5
2026-05-13 14:22:11.038 INFO  [ph2run] *** Circuit: trott ***
2026-05-13 14:22:11.038 INFO  [ph2run] running simulation: 4 Trotter steps
2026-05-13 14:22:11.041 INFO  [trott]  step 1/4 (25%) elapsed 0.00s eta 0.01s
2026-05-13 14:22:11.044 INFO  [trott]  step 2/4 (50%) elapsed 0.01s eta 0.01s
2026-05-13 14:22:11.047 INFO  [trott]  step 3/4 (75%) elapsed 0.01s eta 0.00s
2026-05-13 14:22:11.050 INFO  [trott]  step 4/4 (100%) elapsed 0.01s eta 0.00s
2026-05-13 14:22:11.051 INFO  [ph2run] simulation finished in 0.013 s
```

**Debug capture (make debug build, PHASE2_LOG=debug).**

```sh
make debug
PHASE2_LOG=debug ./build/ph2run/ph2run -S test/data/H2O_CAS56.h5 trott -s 4
```

Adds engine-layer DEBUG lines: `data_open`, qubit layout,
hamil sort, per-step debug markers, cache flush counts.

**Multi-rank capture (PHASE2_LOG_ALL).**

```sh
mpirun -n 4 -x PHASE2_LOG_ALL=1 -x PHASE2_LOG=info \
    ./build/ph2run/ph2run -S test/data/H2O_CAS56.h5 trott -s 4
```

Rank-0 lines remain unprefixed; ranks 1-3 prefix every
line with `[r=N] `.  Pipe to `sort -k1` for chronological
interleaving by rank.
