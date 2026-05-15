# phase2

Full-state vector quantum Hamiltonian simulation library for
quantum phase estimation.

Version 1.1.0

Copyright (c) 2025, Marek Miller.  BSD 3-Clause License.

---

## 1. Introduction

phase2 simulates the time-evolution operator exp(-iHt) for
Hamiltonians H expressed as weighted sums of Pauli strings:

    H = sum_k  c_k  P_k

where each P_k is a tensor product of single-qubit Pauli
operators (I, X, Y, Z) and c_k is a real coefficient.  The
library performs full-state vector simulation on classical
hardware, targeting quantum phase estimation (QPE) workflows
in computational chemistry and quantum algorithm research.

The primary use case is computing the overlap

    <psi| exp(-i H delta)^s |psi>

for a reference state |psi> (given as a multideterminant
expansion) across multiple Trotter steps s, yielding
complex-valued time-series data from which energy estimates
are extracted via classical post-processing.

**Target audience.** Computational chemists studying
molecular electronic structure via QPE; quantum algorithm
researchers benchmarking product formulae and randomised
channel methods on classically tractable system sizes.

**Supported platforms.**

- **CPU**: any POSIX system with MPI.  Tested on GNU/Linux
  with OpenMPI.  The number of MPI ranks must be a power of
  two.
- **GPU**: NVIDIA GPUs via CUDA, combined with MPI for
  multi-node execution.  Requires CUDA-aware MPI (e.g.
  OpenMPI built with CUDA support).

---

## 2. Quick Start

### 2.1 Dependencies (CPU build)

Required packages on Ubuntu/Debian:

    sudo apt install gcc make libopenmpi-dev       \
                     libhdf5-dev

### 2.2 Build and test

    git clone <repository-url> phase2
    cd phase2
    make build
    make build-test
    make check

To run the test suite under MPI with 2 ranks:

    make check-mpi MPIRANKS=2

### 2.3 CUDA build

Additional dependencies:

    sudo apt install nvidia-cuda-toolkit

Or install the CUDA toolkit from NVIDIA directly.  Ensure
`nvcc` is on `PATH` and that the MPI installation is
CUDA-aware.

Build with the CUDA backend:

    make build BACKEND=cuda

The Makefile variable `CUDA_PREFIX` defaults to
`/usr/local/cuda`.  Override if the toolkit is installed
elsewhere:

    make build BACKEND=cuda CUDA_PREFIX=/opt/cuda-12.6

---

## 3. Rationale and Design

### 3.1 Pauli String Encoding

Each Pauli string over up to 64 qubits is encoded in
`struct paulis`, which contains two `uint64_t` fields
`pak[0]` and `pak[1]`.  For qubit n, the single-qubit
Pauli operator is encoded as:

| Operator | pak[0] bit n | pak[1] bit n |
|----------|-------------|-------------|
| I        | 0           | 0           |
| X        | 1           | 0           |
| Z        | 0           | 1           |
| Y        | 1           | 1           |

This is deliberately NOT the natural I=0, X=1, Y=2, Z=3
ordering.  The encoding is I=00, X=10, Z=01, Y=11.  The
rationale:

- `pak[0]` is the "flip mask": bit n is set if and only if
  the operator on qubit n flips the computational basis
  state (X and Y flip; I and Z do not).
- `pak[1]` is the "phase mask": bit n is set if and only if
  the operator on qubit n produces a sign or phase when
  acting on |1> (Z and Y produce -1 and -i respectively;
  I and X do not).

This decomposition enables the core `paulis_effect`
function to be implemented entirely with bitwise operations
and popcount.

**`paulis_effect(code, i, &z)`**: Given a Pauli string P
(encoded as `code`) and a computational basis state |i>,
computes j and a complex phase z such that P|i> = z|j>.

The computation proceeds as:

1. **Bit flip**: j = i XOR pak[0].  The X and Y components
   flip their respective qubits.

2. **Phase accumulation**: The phase is a 4th root of unity
   determined by two popcount operations:
   - mi = popcount(i AND pak[1]) counts the number of
     qubits where |1> is acted on by Z or Y.  Each such
     qubit contributes a factor of -1 (from Z) or a factor
     whose imaginary part contributes (from Y).
   - is = popcount(pak[0] AND pak[1]) counts the number of
     Y operators in the string.  Each Y contributes a
     factor of i (the imaginary unit).
   - The combined phase index is (is + 2*mi) mod 4,
     selecting from {1, i, -1, -i}.

The function multiplies the caller-supplied `*z` by the
computed phase and returns j.  If `z` is NULL, only the
bit-flip index j is returned (used when only the partner
index is needed).

### 3.2 MPI Parallelism Model

The quantum register of N qubits has 2^N amplitudes.  These
are distributed across MPI ranks as follows:

- N = qb_lo + qb_hi, where 2^qb_hi equals the number of
  MPI ranks.
- Each rank stores 2^qb_lo complex double amplitudes in a
  contiguous array `amp[]`.
- Rank r holds the amplitudes for basis states
  |i> where the upper qb_hi bits of i equal r.

**Constraint**: the number of MPI ranks MUST be a power of
two.  The library computes qb_hi = log2(size) during
`qreg_init` and fails if size is not a power of two.

### 3.3 Pauli Rotation Protocol

The core computational kernel applies the unitary
exp(i*phi*P) to the state vector, where P is a Pauli
string and phi is a real angle.  The Pauli string P is
split into hi-qubit and lo-qubit parts:

    P = P_hi (x) P_lo

where P_hi acts on the qb_hi distributed qubits and P_lo
acts on the qb_lo local qubits.

The protocol has four stages:

**1. Exchange.**  P_hi determines which remote rank holds
the partner amplitudes.  `paulis_effect` applied to the
rank index with the hi-part code gives the partner rank.
Non-blocking MPI `Isend`/`Irecv` exchanges the full local
amplitude array.  If the local array exceeds MAX_COUNT
(2^29) elements, the transfer is split into multiple
non-blocking message pairs.

**2. Mix.**  After the exchange completes (MPI_Waitall),
`amp[i]` holds local amplitudes and `buf[i]` holds the
partner rank's amplitudes.  The hi-part phase bm is
computed via `paulis_effect(code_hi, partner_rank, &bm)`.
The mix kernel transforms each element:

    buf[i] *= bm
    x = amp[i],  y = buf[i]
    amp[i] = (x + y) / 2
    buf[i] = (x - y) / 2

This is the exchange-and-split step that separates the
amplitudes into two halves that will rotate in opposite
directions.

**3. Rotate.**  For each pair (i, j) where
P_lo|i> = z|j> with j >= i (to avoid double-counting):

    amp[i] = cos(phi)*amp[i] + i*conj(z)*sin(phi)*amp[j]
    amp[j] = cos(phi)*amp[j] + i*z*sin(phi)*amp[i]

This 2x2 rotation is applied to `amp[]` with angle +phi
and to `buf[]` with angle -phi.

**4. Add.**  The two halves are recombined:

    amp[i] += buf[i]

### 3.4 Cache Batching (circ_cache)

Consecutive Hamiltonian terms that share the same hi-qubit
Pauli string require the same MPI exchange pattern.  Since
MPI exchange is the dominant cost, the library batches such
terms to amortise the communication overhead.

The `circ_cache` module is a static accumulator that stores
up to CACHE_MAX = 1024 (lo-qubit code, angle) pairs.
Terms are inserted via `circ_cache_insert`.  When a term
with a different hi-part arrives, or the cache is full, the
accumulated batch is flushed: a single MPI exchange is
performed, followed by CACHE_MAX (or fewer) local rotation
passes.

Sorting the Hamiltonian lexicographically
(`circ_hamil_sort_lex`) before simulation maximises the
number of consecutive terms sharing the same hi-part,
thereby maximising cache hits and minimising MPI exchanges.

### 3.5 phase2 vs ph2run -- layering and dependencies

phase2 is a compute library: in-memory data carriers, MPI
register, Pauli rotation kernels, and the step-callback
interface that the algorithms (`circ/{trott,trott2,qdrift,
cmpsit}.c`) consume.  It does not link HDF5.

ph2run is the CLI driver: argument parsing, file I/O,
per-step writing to `simul.h5`, summary CSV.  It links
HDF5 and depends on phase2.  No phase2 source file
includes `ph2run/data.h`; the dependency goes one way:

```
+--------- libphase2.so (no HDF5) ---------+
| phase2/circ.c                            |
| phase2/state_prep_coeff.c                |
| phase2/qreg.c / paulis.c / prob.c / ...  |
| circ/{trott, trott2, qdrift, cmpsit}.c   |
| -- phase2/step_writer.h (callback ABI)   |
+------------------------------------------+
                ^
                |  uses
                |
+-- ph2run (CLI; links HDF5 and libphase2.so) --+
| ph2run/ph2run.c                                 |
| ph2run/data.c    -- simul.h5 reads + writes     |
| -- ph2run/data.h: data_open/close,              |
|        data_attr_*, data_state_prep_kind,       |
|        circ_hamil_load, circ_muldet_load,       |
|        data_coeff_matrix_load, data_circ_writer |
+-------------------------------------------------+
```

The carrier types (`struct circ_hamil`, `struct
circ_muldet`, `struct data_coeff_matrix`) live in
`include/phase2/circ.h` because they are the in-memory
shapes the compute kernels consume.  ph2run's loaders
populate them directly: `circ_hamil_load` reads the
on-disk byte arrays from rank 0, broadcasts them, then
packs each row into the `struct circ_hamil_term` form
(coefficient times `norm`, packed `struct paulis`) that
`circ_step` operates on.

### 3.6 Step-writer callback

Per-step output goes through a single function pointer:

```c
/* include/phase2/step_writer.h */
struct phase2_step_writer {
        void *ctx;
        int (*write)(void *ctx, size_t step_idx,
                _Complex double z);
};
```

Each algorithm's `*_init` accepts a `struct
phase2_step_writer *sw`.  Inside `*_simul`, after
computing each step's overlap, the algorithm calls
`sw->write(sw->ctx, i, vals->z[i])`.  Passing `NULL`
disables per-step output and runs the simulation entirely
in memory.

ph2run plugs in a writer whose `ctx` is a
`struct data_circ_writer` (the HDF5 dataset cache from
the data subsystem) and whose `write` thunks to
`data_circ_write_step`:

```c
static int step_thunk(void *ctx, size_t i,
        _Complex double z)
{
        return data_circ_write_step(ctx, i, z);
}
struct phase2_step_writer sw = {
        .ctx = &io.wr, .write = step_thunk };
```

External wrappers (Python via ctypes, future C/Rust
callers) plug in whatever sink they want, without
involving HDF5 or `simul.h5`.

### 3.7 Backend Abstraction

The compile-time macro `PHASE2_BACKEND` selects the
computational backend:

| Value | Macro             | Backend        | Source files              |
|-------|-------------------|----------------|---------------------------|
| 0     | `WORLD_BACKEND "qreg"` | CPU     | qreg_qreg.c              |
| 2     | `WORLD_BACKEND "CUDA"` | NVIDIA GPU | qreg_cuda.c, qreg_cuda_lo.cu |

Both backends implement the same internal interface:

- `qreg_backend_init(reg)`: allocate backend-specific
  resources (e.g. GPU device memory).
- `qreg_backend_free(reg)`: release backend resources.
- `qreg_getamp(reg, i, z)`: retrieve amplitude at global
  index i.
- `qreg_setamp(reg, i, z)`: set amplitude at global
  index i.
- `qreg_zero(reg)`: zero the entire register.
- `qreg_paulirot(reg, code_hi, codes_lo, phis, ncodes)`:
  apply a batch of Pauli rotations sharing the same
  hi-part.

The backend is selected at build time via `make
BACKEND=qreg` (default) or `make BACKEND=cuda`.

---

## 4. Computational Kernels

### 4.1 CPU Kernels (qreg_qreg.c)

The CPU backend implements three inline kernels, each
operating on a single amplitude index i:

**`kernel_mix(i, a, b, bm)`**: Multiplies `b[i]` by the
hi-part phase `bm`, then splits `a[i]` and `b[i]` into
sum and difference halves:

    b[i] *= bm
    x = a[i], y = b[i]
    a[i] = (x + y) / 2
    b[i] = (x - y) / 2

**`kernel_rot(i, a, code, c, s)`**: Computes
`paulis_effect(code, i, &sz)` to find the partner index j
and phase sz.  If j < i, returns immediately (the pair is
handled when the loop reaches the smaller index).
Otherwise applies the 2x2 rotation:

    a[i] = c*a[i] + i*conj(sz)*a[j]
    a[j] = c*a[j] + i*sz*a[i]

where c = cos(phi) and s = sin(phi), with sz initialised
to s before `paulis_effect` modifies it.

**`kernel_add(i, a, b)`**: Accumulates partner amplitudes
back: `a[i] += b[i]`.

The loop structure is simple: for each kernel, a `for` loop
iterates over all `namp` local amplitudes.  For rotation
batches, an outer loop over `ncodes` wraps the inner
amplitude loop.

### 4.2 CUDA Kernels (qreg_cuda_lo.cu)

The CUDA backend maps each amplitude to one GPU thread.
All kernels use 512 threads per block, with
grid = ceil(namp / 512).

**`kernelMix`**: Identical logic to the CPU `kernel_mix`,
using `cuCmul`, `cuCadd`, `cuCsub` for complex arithmetic.

**`paulisEffect` (__device__)**: Device-side implementation
of `paulis_effect`, using `__popcll` (CUDA 64-bit
popcount intrinsic) instead of `stdc_count_ones_ul`.

**`kernelPauliRot`**: One thread per amplitude.  The j < i
guard prevents double-application.  Uses `cuCmul`,
`cuCadd`, `cuConj` for complex operations.

**`kernelAdd`**: One thread per amplitude, `a[i] += b[i]`.

The host function `qreg_paulirot_lo` launches kernels
sequentially in the default CUDA stream.  Kernels within
the same stream execute in order, so no explicit
synchronisation is needed between them.

### 4.3 MPI Communication (qreg.c, qreg_cuda.c)

**Message chunking.**  The maximum number of complex double
elements per MPI message is MAX_COUNT = 2^29.  If the
local amplitude array exceeds this limit, it is
partitioned into `nreqs = namp / msg_count` chunks, each
transferred via a separate `MPI_Isend`/`MPI_Irecv` pair.

**CUDA-aware MPI.**  The CUDA backend passes device
pointers directly to `MPI_Isend` and `MPI_Irecv`.  This
requires a CUDA-aware MPI implementation that can
transfer data directly from GPU memory via RDMA or staging
buffers.  A `cudaDeviceSynchronize` call precedes the MPI
operations to ensure all pending GPU work has completed.

---

## 5. Algorithms

### 5.1 Trotter (circ/trott.c)

First-order Lie-Trotter product formula.  For a
Hamiltonian H = sum_k c_k P_k, each Trotter step applies:

    prod_k exp(i * delta * c_k * P_k)

sequentially for each term k.  The parameter `delta` is
the rotation angle (related to the time step by the
Hamiltonian normalisation).

The Hamiltonian is sorted lexicographically on
initialisation (`circ_hamil_sort_lex`) to maximise
circ_cache batch efficiency.

The simulation loop:

1. Prepare initial state from the multideterminant
   expansion.
2. For each step s = 1, ..., steps:
   a. Apply one Trotter step via `circ_step`.
   b. Measure the overlap <psi|state> and store in
      `vals[s]`.

### 5.2 qDRIFT (circ/qdrift.c)

Randomised product formula based on the qDRIFT channel
(Campbell, 2019).  Instead of applying all Hamiltonian
terms deterministically, each step samples terms from the
Hamiltonian with probability proportional to |c_k|.

**Initialisation.**  A cumulative distribution function
(CDF) is built from the absolute values of the Hamiltonian
coefficients |c_k|, normalised to a probability
distribution.

**Simulation loop.**  For each of `samples` independent
samples:

1. Prepare the initial state.
2. Draw `depth` terms independently from the CDF.  For
   each drawn term k, apply:

       exp(i * asin(step_size) * sign(c_k) * P_k)

   The `asin` transformation relates the `step_size`
   parameter to the physical rotation angle.
3. Measure the overlap and store the result.

The PRNG is xoshiro256** (Blackman and Vigna, 2021),
initialised from the user-supplied seed.  The PRNG state
is deterministically split across MPI ranks via the
`world_init` seed-splitting mechanism.

### 5.3 Composite (circ/cmpsit.c)

Partially randomised second-order Suzuki-Trotter formula.
The Hamiltonian is split into two parts:

- **Deterministic part**: the top L terms by |c_k|
  (parameter `length`), sorted lexicographically for cache
  efficiency.
- **Randomised part**: the remaining terms, sampled via
  CDF (same mechanism as qDRIFT).

Each Trotter step consists of:

1. Sample a composite Hamiltonian: L deterministic terms
   (with coefficients scaled by `angle_det`) followed by
   `depth` randomly sampled terms (with rotation angle
   `angle_rand` times sign(c_k)).
2. Forward half-step: `circ_step` with omega = 0.5.
3. Fresh sample (independent draw).
4. Reverse half-step: `circ_step_reverse` with omega = 0.5
   (terms applied in reverse order, completing the
   second-order symmetric decomposition).

The reverse-order application in step 4 is the key
difference from a naive first-order scheme: it yields the
symmetric S2 integrator, which has second-order error
scaling.

---

## 6. Usage

### 6.1 The `ph2run` CLI

    ph2run [OPTIONS] CMD [CMD_OPTIONS]

**Global options:**

| Flag        | Description                              |
|-------------|------------------------------------------|
| `-v`        | Verbose output                           |
| `-S FILE`   | Simulation HDF5 file (default: ./simul.h5) |
| `--version` | Print version and exit                   |
| `--help`    | Print help and exit                      |

### 6.2 Subcommands

#### `trott` -- 1st-order Trotter

    ph2run [OPTS] trott [TROTT_OPTS]

| Flag    | Description                          | Default |
|---------|--------------------------------------|---------|
| `-D VAL`| Trotter step size (delta)            | 1.0     |
| `-s N`  | Number of Trotter steps              | 1       |

#### `trott2` -- 2nd-order symmetric (Strang) Trotter

    ph2run [OPTS] trott2 [TROTT2_OPTS]

| Flag    | Description                          | Default |
|---------|--------------------------------------|---------|
| `-D VAL`| Trotter step size (delta)            | 1.0     |
| `-s N`  | Number of Trotter steps              | 1       |

#### `qdrift` -- qDRIFT randomised product formula

    ph2run [OPTS] qdrift [QDRIFT_OPTS]

| Flag    | Description                          | Default |
|---------|--------------------------------------|---------|
| `-D VAL`| qDRIFT step size                     | 1.0     |
| `-d N`  | Randomised terms per sample          | 64      |
| `-n N`  | Number of independent samples        | 1       |
| `-x N`  | PRNG seed (must be non-zero)         | 1       |

#### `cmpsit` -- Composite (2nd-order, partially randomised)

    ph2run [OPTS] cmpsit [CMPSIT_OPTS]

| Flag    | Description                          | Default |
|---------|--------------------------------------|---------|
| `-l N`  | Number of deterministic top-|c_k|    | 1       |
| `-d N`  | Randomised terms per step            | 64      |
| `-s N`  | Number of Trotter steps              | 1       |
| `-D VAL`| Deterministic step size (angle_det)  | 1.0     |
| `-R VAL`| Randomised step size (angle_rand)    | 1.0     |
| `-n N`  | Number of independent samples        | 1       |
| `-x N`  | PRNG seed (must be non-zero)         | 1       |

### 6.3 Environment Variables

- **`PHASE2_LOG`** (`info` default): one of `trace`,
  `debug`, `info`, `warn`, `error`, `fatal`.  Raises
  verbosity above the default.
- **`PHASE2_LOG_ALL`** (unset default): any non-empty
  value opens logging on rank > 0 in MPI runs; rank > 0
  lines are prefixed `[r=N]`.

Defaults: `info`-and-below on stdout, `warn`-and-above on
stderr, both line-buffered + flushed per emit (SLURM
`tail -f` friendly).  Release builds strip `log_trace`
and `log_debug` to no-ops — rebuild with `make debug` to
trace a production-class run.

See [doc/logging.md](logging.md) for the full reference
(format, policy for new code, performance budget,
multi-rank semantics, examples).

---

## 7. Examples

The test suite includes several HDF5 input files:

- `test/data/H2O_CAS56.h5`: water in a CAS(5,6) active space,
  10 qubits, 251 Hamiltonian terms, single-determinant
  reference (`/state_prep/multidet`).
- `test/data/case-d9f603dc.h5_solved`: 3-qubit toy case,
  10 terms, 3 determinants.
- `test/data/N4_closed.h5`: 8-qubit, n_sites=4, closed-shell
  reference encoded as `/state_prep/coeff_matrix` — exercises
  the Slater-Condon expansion path.
- `test/data/bendazzoli/n4_oss_k0/n4_oss_k0.h5`: open-shell
  CSF superposition (two-component) used by the
  `t-ref-bendazzoli` reference test.

### 7.1 Trotter (4 steps, 2 MPI ranks)

    mpirun -n 2 ./ph2run/ph2run -S test/data/H2O_CAS56.h5 \
        trott -D 0.1 -s 4

### 7.2 Symmetric (Strang) Trotter

    mpirun -n 2 ./ph2run/ph2run -S test/data/H2O_CAS56.h5 \
        trott2 -D 0.1 -s 4

### 7.3 qDRIFT (depth 64, 100 samples)

    mpirun -n 2 ./ph2run/ph2run -S test/data/H2O_CAS56.h5 \
        qdrift -D 0.05 -d 64 -n 100 -x 42

### 7.4 Composite (2nd-order)

    mpirun -n 2 ./ph2run/ph2run                            \
        -S test/data/case-d9f603dc.h5_solved               \
        cmpsit -l 3 -d 32 -s 4 -D 0.1 -R 0.05 -n 50 -x 7

### 7.5 Coefficient-matrix state prep

    mpirun -n 2 ./ph2run/ph2run -S test/data/N4_closed.h5 \
        trott -D 0.05 -s 8

The reference state is reconstructed at `circ_prepst` time
from `/state_prep/coeff_matrix`; no per-determinant
amplitudes appear in the file (see §3 of
`doc/simul-h5-specs.md`).

All commands write results back into the respective HDF5
groups in the simulation file.

---

## 8. API Reference

### 8.1 Pauli Strings (`include/phase2/paulis.h`)

```c
struct paulis {
    uint64_t pak[2];
};
```

Packed representation of a Pauli string over up to 64
qubits.  See Section 3.1 for encoding details.

---

```c
struct paulis paulis_new(void);
```

Return a new Pauli string initialised to the identity on
all qubits (pak[0] = pak[1] = 0).

---

```c
int paulis_get(struct paulis code, uint32_t n);
```

Return the single-qubit Pauli operator at qubit n.

**Parameters:**
- `code`: the Pauli string.
- `n`: qubit index, 0 <= n < 64.

**Return value:** One of PAULI_I (0), PAULI_X (1),
PAULI_Y (2), PAULI_Z (3).

---

```c
void paulis_set(struct paulis *code, int op, uint32_t n);
```

Set the single-qubit Pauli operator at qubit n.

**Parameters:**
- `code`: pointer to the Pauli string to modify.
- `op`: one of PAULI_I, PAULI_X, PAULI_Y, PAULI_Z.
- `n`: qubit index, 0 <= n < 64.

---

```c
int paulis_eq(struct paulis code1, struct paulis code2);
```

Test two Pauli strings for equality.

**Return value:** Nonzero if equal, zero otherwise.

---

```c
void paulis_shl(struct paulis *code, uint32_t n);
```

Left-shift the Pauli string by n qubit positions.
Equivalent to tensoring n identity operators on the right.

---

```c
void paulis_shr(struct paulis *code, uint32_t n);
```

Right-shift the Pauli string by n qubit positions.
Discards the lowest n qubit operators.

---

```c
uint64_t paulis_effect(struct paulis code, uint64_t i,
                       _Complex double *z);
```

Compute the action of Pauli string P on basis state |i>.

**Parameters:**
- `code`: the Pauli string P.
- `i`: the basis state index.
- `z`: pointer to a complex number.  On entry, *z holds
  an input factor.  On exit, *z is multiplied by the phase
  of P|i>.  May be NULL, in which case only the output
  index is computed.

**Return value:** The index j such that P|i> = z|j>.

**Preconditions:** i must be a valid basis state index for
the number of qubits encoded in code.

---

```c
void paulis_split(struct paulis code, uint32_t qb_lo,
    uint32_t qb_hi, struct paulis *lo, struct paulis *hi);
```

Split a Pauli string into lo-qubit and hi-qubit parts.

**Parameters:**
- `code`: the Pauli string to split.
- `qb_lo`: number of local (low) qubits.
- `qb_hi`: number of distributed (high) qubits.
- `lo`: output, receives the lo-qubit part.
- `hi`: output, receives the hi-qubit part (bits remain
  in their original positions, not shifted).

---

```c
void paulis_merge(struct paulis *code, uint32_t qb_lo,
    uint32_t qb_hi, struct paulis lo, struct paulis hi);
```

Merge lo-qubit and hi-qubit parts into a single Pauli
string.  Inverse of `paulis_split`.

---

```c
int paulis_cmp(struct paulis a, struct paulis b);
```

Lexicographic comparison of two Pauli strings, comparing
from the highest qubit index downward.

**Return value:** -1 if a < b, 0 if a == b, 1 if a > b.

---

### 8.2 Quantum Register (`include/phase2/qreg.h`)

```c
#define QREG_MAX_WIDTH (64)

struct qreg {
    struct world_info wd;
    uint32_t qb_lo, qb_hi;
    _Complex double *amp, *buf;
    uint64_t namp;
    int msg_count;
    MPI_Request *reqs_snd, *reqs_rcv;
    size_t nreqs;
    void *data;
};
```

Distributed quantum register.  `amp` points to the local
amplitude array of size `namp = 2^qb_lo`.  `buf` is a
communication buffer of the same size, allocated
contiguously after `amp`.  `data` is an opaque handle to
backend-specific resources (e.g. GPU device pointers for
the CUDA backend).

---

```c
int qreg_init(struct qreg *reg, uint32_t qb);
```

Initialise a distributed quantum register of `qb` qubits.

**Parameters:**
- `reg`: pointer to an uninitialised qreg struct.
- `qb`: total number of qubits.

**Return value:** 0 on success, -1 on error.

**Preconditions:**
- `world_init` must have been called.
- The number of MPI ranks must be a power of two.
- qb must be greater than qb_hi = log2(ranks).

---

```c
void qreg_free(struct qreg *reg);
```

Release all resources associated with the quantum
register, including backend resources, amplitude arrays,
and MPI request arrays.

---

```c
void qreg_getamp(struct qreg *reg, uint64_t i,
                 _Complex double *z);
```

Retrieve the amplitude at global basis state index i.

**Parameters:**
- `reg`: initialised quantum register.
- `i`: global basis state index, 0 <= i < 2^(qb_lo+qb_hi).
- `z`: output, receives the amplitude value.

This is a collective MPI operation (broadcasts the value
from the owning rank to all ranks).

---

```c
void qreg_setamp(struct qreg *reg, uint64_t i,
                 _Complex double z);
```

Set the amplitude at global basis state index i.

**Parameters:**
- `reg`: initialised quantum register.
- `i`: global basis state index.
- `z`: the amplitude value to set.

This is a collective MPI operation (barrier after the
set).

---

```c
void qreg_zero(struct qreg *reg);
```

Set all amplitudes to zero.  Collective MPI operation.

---

```c
void qreg_paulirot(struct qreg *reg,
    struct paulis code_hi,
    const struct paulis *codes_lo,
    const double *phis, size_t ncodes);
```

Apply a batch of Pauli rotations sharing the same hi-qubit
part.  For each k in [0, ncodes), applies
exp(i * phis[k] * (code_hi (x) codes_lo[k])) to the
register.

**Parameters:**
- `reg`: initialised quantum register.
- `code_hi`: the shared hi-qubit Pauli string (bits in
  original positions, not shifted).
- `codes_lo`: array of ncodes lo-qubit Pauli strings.
- `phis`: array of ncodes rotation angles.
- `ncodes`: number of rotations in the batch.

**Preconditions:** All codes_lo entries must have zero bits
in the hi-qubit positions.

---

### 8.3 World (`include/phase2/world.h`)

```c
enum world_stat {
    WORLD_UNDEF = -1,
    WORLD_READY = 0,
    WORLD_DONE  = 1,
    WORLD_ERR   = 2,
};

struct world_info {
    enum world_stat stat;
    int size;
    int rank;
    uint64_t seed;
};
```

---

```c
int world_init(int *argc, char ***argv, uint64_t seed);
```

Initialise the global simulation environment.  Must be
called exactly once at the start of the program.
Initialises MPI, the logging facility, and the PRNG.

**Parameters:**
- `argc`, `argv`: pointers to command-line argument count
  and vector (forwarded to MPI_Init).  May be NULL.
- `seed`: seed for the xoshiro256** PRNG.  Must not be
  zero.  The PRNG state is deterministically split across
  MPI ranks.

**Return value:** WORLD_READY on success, WORLD_ERR on
failure.

---

```c
int world_free(void);
```

Destroy the global simulation environment.  Finalises MPI
and the logging facility.  Must be called exactly once at
the end of the program.

**Return value:** WORLD_DONE on success, WORLD_ERR on
failure.

---

```c
int world_info(struct world_info *wd);
```

Populate `wd` with information about the current world
state (rank, size, seed, status).

**Return value:** The value of `wd->stat` after the call.

---

### 8.4 Data I/O (`include/ph2run/data.h`)

The data subsystem lives on the ph2run side and is the
only HDF5-aware piece of the build.  phase2 sources do
not include this header.  The carrier types it populates
(`struct circ_hamil`, `struct circ_muldet`, `struct
data_coeff_matrix`) are defined in
`include/phase2/circ.h`.

```c
typedef int64_t data_id;
#define DATA_INVALID_FID   INT64_C(-1)
#define DATA_FOLLOWER_FID  INT64_C(-2)
```

`DATA_FOLLOWER_FID` is the rank > 0 sentinel returned by
`data_open`: the data layer is rank-0-only; non-zero
ranks short-circuit reads (receiving the bcast result)
and writes (no-ops).

---

```c
data_id data_open(const char *filename);
void    data_close(data_id fid);
```

`data_open` returns a positive file id on rank 0,
`DATA_FOLLOWER_FID` on follower ranks, or
`DATA_INVALID_FID` on error.  Both calls are collective.

---

```c
int data_grp_create(data_id fid, const char *grp_name);
```

Idempotently create an HDF5 group at the root of the
file.  Handles dangling soft-links left by older
fixtures.

---

```c
int  data_attr_read_dbl(data_id fid, const char *grp,
        const char *name, double *a);
int  data_attr_write_dbl(data_id fid, const char *grp,
        const char *name, double a);
int  data_attr_write_ul(data_id fid, const char *grp,
        const char *name, unsigned long a);

#define data_attr_read(fid, grp, name, ptr)  /* _Generic dispatch */
#define data_attr_write(fid, grp, name, val) /* _Generic dispatch */
```

The `_Generic` dispatchers route by argument type.  Read
covers `double *`; write covers `double` and `unsigned
long`.

---

```c
enum stprep_kind {
        STPREP_MULTIDET     = 1,
        STPREP_COEFF_MATRIX = 2,
};

int data_state_prep_kind(data_id fid,
        enum stprep_kind *out);
```

Probe the file for `/state_prep/multidet` vs
`/state_prep/coeff_matrix`.  Both-present and
neither-present are errors (rebuild simul.h5).  See
`doc/simul-h5-specs.md` for the dispatch table.

---

```c
int circ_hamil_load(data_id fid, struct circ_hamil *hm);
int circ_muldet_load(data_id fid, struct circ_muldet *md);
int data_coeff_matrix_load(data_id fid,
        struct data_coeff_matrix *cm);
```

Read a `/pauli_hamil`, `/state_prep/multidet`, or
`/state_prep/coeff_matrix` group end-to-end: rank 0
reads, all ranks receive via `MPI_Bcast`, and -- for hamil
and multidet -- each rank packs the raw bytes into the
phase2-side struct (`circ_hamil_term`/`circ_muldet`
entries) the compute kernels consume.  Return 0 on
success, -1 on error (with `log_error`).

The matching `*_free` calls live in
`include/phase2/circ.h` (pure pointer cleanup, no HDF5):
`circ_hamil_free`, `circ_muldet_free`,
`data_coeff_matrix_free`.

---

```c
struct data_circ_writer {
        data_id fid;
        int64_t dset;     /* H5Dopen handle on rank 0 */
        size_t  n_steps;
};

int  data_circ_writer_init(data_id fid,
        const char *grp_name, size_t n_steps,
        struct data_circ_writer *w);
int  data_circ_write_step(struct data_circ_writer *w,
        size_t step_idx, _Complex double z);
void data_circ_writer_close(struct data_circ_writer *w);
```

Pre-allocate a NaN-padded `(n_steps, 2)` `values`
dataset under `grp_name` and cache the open dataset
handle in `*w`.  Per-step writes hyperslab one row and
flush so a crash leaves the file consistent up to the
last flushed row.  Passing `fid == 0` to
`data_circ_writer_init` makes subsequent writes no-ops
(for callers that want to run without per-step output).

ph2run wraps this writer in a `struct phase2_step_writer`
(see Section 3.6) and hands it to the algorithm's
`*_init`.

---

### 8.5 Circuit Infrastructure (`include/phase2/circ.h`)

```c
struct circ_hamil {
    uint32_t qb;
    struct circ_hamil_term {
        double cf;
        struct paulis op;
    } *terms;
    size_t len;
};
```

In-memory representation of a Pauli Hamiltonian.

---

```c
struct circ_muldet {
    struct { uint64_t idx; _Complex double cf; } *dets;
    size_t len;
};
```

Multideterminant initial state: a sum of computational
basis states with complex coefficients.

---

```c
struct circ {
    struct circ_hamil hm;
    struct circ_muldet md;
    struct circ_cache *cache;
    struct circ_values vals;
    struct qreg reg;
    enum stprep_kind stprep_kind;
    struct data_coeff_matrix cm;
};
```

Complete simulation context: Hamiltonian, multideterminant
or coefficient-matrix state-prep (selected by
`stprep_kind`), quantum register, cache, and output
values buffer.  The carrier types `struct
data_coeff_matrix` and `enum stprep_kind` are defined in
this same header.

---

```c
int circ_hamil_init(struct circ_hamil *hm, uint32_t qb,
                    size_t len);
```

Allocate a Hamiltonian with `len` terms over `qb` qubits.

**Return value:** 0 on success, -1 on allocation failure.

---

```c
void circ_hamil_free(struct circ_hamil *hm);
```

Free the terms array.

---

```c
void circ_hamil_sort_lex(struct circ_hamil *hm);
```

Sort Hamiltonian terms in lexicographic order of Pauli
strings (highest qubit first).  This maximises circ_cache
batch efficiency.

---

```c
int circ_muldet_init(struct circ_muldet *md, size_t len);
void circ_muldet_free(struct circ_muldet *md);
```

Allocate / free a multideterminant initial state with
`len` determinants.

---

```c
void circ_prog_init(struct circ_prog *prog, size_t len);
void circ_prog_tick(struct circ_prog *prog);
```

Progress tracker.  `circ_prog_init` sets the total number
of steps.  `circ_prog_tick` increments the counter and
emits a log_info message each time a new percentage point
is reached.

---

```c
int circ_values_init(struct circ_values *vals, size_t len);
void circ_values_free(struct circ_values *vals);
```

Allocate / free an output buffer for `len` complex values.

---

```c
int circ_init(struct circ *ct, struct circ_hamil hm,
        enum stprep_kind sp_kind, const void *sp_data,
        size_t vals_len);
```

Adopt a pre-loaded Hamiltonian and state-prep payload
into a fresh circ context, allocate the quantum register,
initialise the cache, and allocate the output buffer.

The caller (typically ph2run) loads `hm` via
`circ_hamil_load` and the state-prep struct via
`circ_muldet_load` or `data_coeff_matrix_load`, then
passes them in.  `circ_init` takes ownership: the
buffers carried by `hm` and `*sp_data` are owned by
`ct` after a successful call and freed by `circ_free`.
On error the function frees whatever it has adopted, so
the caller's locals must not be freed a second time.

`sp_data` points to a `struct circ_muldet` when
`sp_kind == STPREP_MULTIDET`, or a
`struct data_coeff_matrix` when
`sp_kind == STPREP_COEFF_MATRIX`.

**Return value:** 0 on success, -1 on error.

---

```c
void circ_free(struct circ *ct);
```

Free all resources in the simulation context.

---

```c
int circ_prepst(struct circ *ct);
```

Prepare the initial state: zero the register, then set
amplitudes from the multideterminant expansion.

**Return value:** 0 on success.

---

```c
int circ_step(struct circ *ct,
              const struct circ_hamil *hm, double omega);
```

Apply one forward Trotter step: for each term k (in order),
apply exp(i * omega * c_k * P_k).  Uses circ_cache
internally.

**Return value:** 0 on success, -1 on error.

---

```c
int circ_step_reverse(struct circ *ct,
    const struct circ_hamil *hm, double omega);
```

Apply one reverse Trotter step: same as `circ_step` but
terms are applied in reverse order (last to first).  Used
for the second half of a symmetric Suzuki-Trotter
decomposition.

**Return value:** 0 on success, -1 on error.

---

```c
_Complex double circ_measure(struct circ *ct);
```

Compute the overlap <psi|state> where |psi> is the
multideterminant reference state and |state> is the current
register state.

**Return value:** The complex overlap value.

---

### 8.6 Probability (`include/phase2/prob.h`)

```c
struct prob_cdf {
    double *y;
    size_t len;
};
```

Discrete cumulative distribution function.

---

```c
int prob_cdf_init(struct prob_cdf *cdf, size_t len);
```

Allocate a CDF of length `len`.

**Preconditions:** len > 0.

**Return value:** 0 on success, -1 on allocation failure.

---

```c
void prob_cdf_free(struct prob_cdf *cdf);
```

Free the CDF array.

---

```c
int prob_cdf_from_iter(struct prob_cdf *cdf,
    double (*iter)(void *), void *data);
```

Build the CDF from an iterator.  Calls `iter(data)`
exactly `cdf->len` times to obtain samples.  Takes
absolute values, normalises to a probability distribution,
and computes the cumulative sum.

**Return value:** 0 on success.

---

```c
size_t prob_cdf_inverse(const struct prob_cdf *cdf,
                        double y);
```

Inverse CDF lookup.  Returns the largest index i such that
F(i) <= y, where F is the cumulative distribution function.
Returns 0 if y is below F(0).

---

### 8.7 Trotter (`include/circ/trott.h`)

```c
struct trott_data { double delta; size_t steps; };
struct trott {
        struct circ ct;
        struct trott_data dt;
        struct phase2_step_writer *sw;
};
```

---

```c
int trott_init(struct trott *tt,
        const struct trott_data *dt,
        struct circ_hamil hm, enum stprep_kind sp_kind,
        const void *sp_data,
        struct phase2_step_writer *sw);
void trott_free(struct trott *tt);
int  trott_simul(struct trott *tt);
```

`trott_init` adopts a pre-loaded Hamiltonian and
state-prep payload via `circ_init` (Section 8.5), sorts
the Hamiltonian lexicographically, and stores the
step-writer callback.  `sw == NULL` runs without
per-step output.  `trott_simul` runs `dt.steps` Trotter
steps, populates `ct.vals` with overlaps, and forwards
each row through `sw->write` when non-NULL.

---

### 8.8 Symmetric Trotter (`include/circ/trott2.h`)

Strang 2nd-order product formula.  One step is a forward
sweep at `delta/2` followed by a reverse sweep at
`delta/2`, both expressed through `circ_step` /
`circ_step_reverse`.

```c
struct trott2_data { double delta; size_t steps; };
struct trott2 {
        struct circ ct;
        struct trott2_data dt;
        struct phase2_step_writer *sw;
};

int  trott2_init(struct trott2 *t2,
        const struct trott2_data *dt,
        struct circ_hamil hm, enum stprep_kind sp_kind,
        const void *sp_data,
        struct phase2_step_writer *sw);
void trott2_free(struct trott2 *t2);
int  trott2_simul(struct trott2 *t2);
```

Behaviour mirrors `trott_*`.  See
`include/circ/trott2.h` for the per-function contract.

---

### 8.9 qDRIFT (`include/circ/qdrift.h`)

```c
struct qdrift_data {
        size_t depth;
        size_t samples;
        double step_size;
        uint64_t seed;
};

struct qdrift {
        struct circ ct;
        struct qdrift_data dt;
        struct qdrift_ranct ranct;
        struct xoshiro256ss rng;
        struct phase2_step_writer *sw;
};

int qdrift_init(struct qdrift *qd,
        const struct qdrift_data *dt,
        struct circ_hamil hm, enum stprep_kind sp_kind,
        const void *sp_data,
        struct phase2_step_writer *sw);
void qdrift_free(struct qdrift *qd);
int  qdrift_simul(struct qdrift *qd);
```

`qdrift_init` adopts the carriers, builds the CDF over
the Hamiltonian's |c_k|, and seeds the PRNG.  When
`dt.seed == 0` the function resolves the seed to the
compiled-in default and writes the chosen value back
into `qd->dt.seed`, so the caller can record it on
disk.

---

### 8.10 Composite (`include/circ/cmpsit.h`)

```c
struct cmpsit_data {
        uint64_t seed;
        size_t length;
        size_t depth;
        size_t steps;
        double angle_det;
        double angle_rand;
        size_t samples;
};

struct cmpsit {
        struct circ ct;
        struct cmpsit_data dt;
        struct cmpsit_ranct ranct;
        struct xoshiro256ss rng;
        struct phase2_step_writer *sw;
};

int cmpsit_init(struct cmpsit *cp,
        const struct cmpsit_data *dt,
        struct circ_hamil hm, enum stprep_kind sp_kind,
        const void *sp_data,
        struct phase2_step_writer *sw);
void cmpsit_free(struct cmpsit *cp);
int  cmpsit_simul(struct cmpsit *cp);
```

Same shape as the others.  `cmpsit_init` splits the
Hamiltonian into deterministic and randomised parts,
sorts the deterministic part lexicographically, builds
the CDF for the randomised part, and seeds the PRNG
(with the same `seed == 0` resolution as qDRIFT).
`cmpsit_simul` runs `dt.samples` independent samples,
each of `dt.steps` second-order Trotter steps.

---

## 9. Data Format

All simulation input and output is stored in HDF5 files.
The canonical simulation file is `simul.h5`, whose
structure is fully specified in
[doc/simul-h5-specs.md](simul-h5-specs.md).

Key groups:

- `/state_prep/multidet`: initial state as an explicit
  multideterminant expansion (complex coefficients and
  Slater determinants).
- `/state_prep/coeff_matrix`: initial state as a real
  `(n_sites, n_occ)` molecular-orbital coefficient matrix.
  The simulator expands it to a dense superposition at
  `circ_prepst` time via the Slater-Condon formula.
  Supports `closed_shell`, `tapered`, and an optional
  `csf/` subgroup for CSF superpositions.
- Exactly one of `/state_prep/multidet` or
  `/state_prep/coeff_matrix` must be present; the
  dispatch table is documented in
  [doc/simul-h5-specs.md](simul-h5-specs.md).
- `/pauli_hamil`: Hamiltonian as a list of real
  coefficients and Pauli strings, with a normalisation
  factor.
- `/circ_trott`: 1st-order Trotter results (delta
  attribute and complex-valued output).
- `/circ_trott2`: symmetric (Strang) 2nd-order Trotter
  results (delta attribute and complex-valued output).
- `/circ_qdrift`: qDRIFT results (step_size, depth,
  num_samples, seed attributes and complex-valued output).
- `/circ_cmpsit`: composite results (length, depth,
  angle_det, angle_rand, steps, seed attributes and
  complex-valued output).

Pauli operators in HDF5 datasets use the standard encoding:
I=0, X=1, Y=2, Z=3.  This differs from the internal packed
encoding in `struct paulis` (Section 3.1); the conversion
happens during data loading in `circ_hamil_load`.

---

## 10. Building and Testing

### 10.1 Make Targets

| Target                 | Description                           |
|------------------------|---------------------------------------|
| `all`                  | Build everything (programs, benchmarks, tests) |
| `build`                | Build `ph2run` and dependencies       |
| `build-test`           | Build the test suite + runner         |
| `build-bench`          | Build benchmarks                      |
| `check`                | Run the full test suite (parallel, cargo-style) |
| `check-mpi`            | Run the suite under `mpirun -n MPIRANKS` |
| `check-slow`           | Run TESTS_SLOW (excluded from default `check`) |
| `check-<filter>`       | Run the subset matching `<filter>*` (`check-data`) |
| `check-tests-coverage` | Guard: every `test/t-*.c` is in TESTS / TESTS_SLOW |
| `check-srcs-coverage`  | Guard: every subsys-dir source is in DECLARED_SRCS |
| `bench`                | Run benchmarks (single rank)          |
| `bench-mpi`            | Run benchmarks under MPI              |
| `debug`                | Build with debug flags (-g -Og)       |
| `clean`                | Remove the `build/` tree              |
| `distclean`            | Also remove binaries and `libphase2.so` |
| `format`               | Run clang-format on all sources       |

### 10.2 Backend Selection

    make build BACKEND=qreg     # CPU (default)
    make build BACKEND=cuda     # NVIDIA GPU

The CUDA backend requires `nvcc` and the CUDA toolkit.
The NVCC flags default to `-O3 -dopt=on -arch=native`.
For specific GPU architectures (e.g. NVIDIA H100), override
NVCCFLAGS:

    make build BACKEND=cuda \
        NVCCFLAGS="-O3 -dopt=on -arch=sm_90a"

### 10.3 MPI Configuration

The number of MPI ranks for `check-mpi` and `bench-mpi` is
controlled by `MPIRANKS` (default: 2).  Additional MPI
runtime flags can be passed via `MPIFLAGS`:

    make check-mpi MPIRANKS=4
    make check-mpi MPIRANKS=8 \
        MPIFLAGS="--oversubscribe"

### 10.4 Test Suite Overview

The C test binaries live in `test/` and follow the prefix
convention `t-<area>[_<aspect>]`:

| Area              | Coverage                                |
|-------------------|-----------------------------------------|
| `t-paulis`        | Pauli string encoding and arithmetic    |
| `t-qreg`, `t-bitstring_index` | quantum register layout, MPI ownership |
| `t-world`         | global init / free lifecycle            |
| `t-data_*`        | HDF5 attr, hamil, multidet, coeff_matrix, trott-step I/O |
| `t-prob`          | CDF construction and inversion          |
| `t-combinations`, `t-det_small` | enumerator and determinant primitives |
| `t-circ`, `t-circ_cache`        | circuit infrastructure and cache batching |
| `t-circ_trott`, `t-circ_trott2` | end-to-end Trotter / Strang Trotter |
| `t-circ_trott_coeff`, `t-circ_trott2_coeff` | same, on coeff_matrix inputs |
| `t-circ_prepst_coeff`           | Slater-Condon state-prep dispatch |
| `t-state_prep_coeff_*`          | expand, CSF superposition, large reference |
| `t-ref-bendazzoli`              | precomputed CSF amplitudes reference |
| `t-run`                         | self-tests for the parallel runner   |

The Python harness `test/t-ref-coeff_matrix.py` cross-checks
the C expansion against an in-tree reference oracle in
`test/ref/coeff_matrix_reference.py`; the unified `check`
target runs it alongside the C tests.

`make check` and friends delegate to the cargo-style
parallel runner at `test/run` -- per-test stdout/stderr
capture, TTY-aware colour, fnmatch filtering, MPI wrap.
See [doc/testing.md](testing.md) for the full subsystem
reference: macro contract, fixture pattern, runner flags,
recipes for adding tests, sanitiser / valgrind targets.
Test data files reside in `test/data/`.

### 10.5 Build Organisation

The build is split into a top-level orchestrator
plus a per-subsystem Makefile in each source dir:

    Makefile                <- orchestrator + config
    build-rules.mk          <- shared compile rules
    phase2/Makefile         <- compiles phase2 sources
    circ/Makefile
    lib/Makefile
    ph2run/Makefile         <- + links ph2run
    bench/Makefile          <- + links b-paulis, b-qreg
    test/Makefile           <- + tests, runner, fixtures

The top-level Makefile reads in three sections.
Section 1 is **USER CONFIGURATION**: toolchain
(`CC`, `CFLAGS`, `EXTRA_CFLAGS`), MPI library paths
(`MPI_CFLAGS`, `MPI_LDFLAGS`, `MPI_LDLIBS`), HDF5
paths (`HDF5_CFLAGS`, `HDF5_LDFLAGS`,
`HDF5_LDLIBS`), backend selection (`BACKEND`,
`CUDA_PREFIX`, `NVCCFLAGS`), `BUILDDIR`, and
`VERSION`.  A user porting phase2 to a new machine
edits this block only.  Section 2 is **INTERNAL**:
subsystem layout, derived OBJS unions, the CFLAGS
pipeline, the `export` contract that makes Section 1
values visible to sub-makes.  Section 3 is
**TARGETS**: per-subsys dispatch (`$(MAKE) -C
phase2`, ...) and the user-visible high-level
targets.

Every compile artefact -- objects, header-dep files,
final binaries, the shared library -- lands under
`$(BUILDDIR)/<srcdir>/`, mirroring the source tree:

    phase2/circ.c    ->  build/phase2/circ.o + build/phase2/circ.d
    lib/log.c        ->  build/lib/log.o     + build/lib/log.d
    ph2run/ph2run.c  ->  build/ph2run/ph2run (binary)
    test/t-paulis.c  ->  build/test/t-paulis (binary)
    bench/b-paulis.c ->  build/bench/b-paulis (binary)
    test/run.c       ->  build/test/run (binary)
    (shared lib)     ->  build/libphase2.so

Source directories stay clean -- no `.o`, no
binaries.  `make clean` is a single `rm -rf build/`.

Header dependencies are tracked automatically via
`gcc -MMD -MP`: every compile emits a
`build/<dir>/<x>.d` describing which headers the
`.c` `#include`'d transitively, with empty fallback
rules for each header so a renamed header doesn't
break the build.  The top-level Makefile
`-include`s the dep files at its foot; the wildcard
is silently empty on first build.  Touching a
header reliably rebuilds every dependent object --
there are no hand-declared per-object header
prereqs to keep in sync.

The shared library (`libphase2.so`) compiles the
same sources with `-fPIC` into a disjoint tree at
`build/shared/<srcdir>/`.  Running `make build` and
`make shared` (in either order) cannot mix PIC and
non-PIC objects, because the two layouts live on
different paths.

Each subsys Makefile carries its own
`check-srcs-coverage` rule: it enumerates `*.c` (and
`*.cu` for phase2) in its own dir and fails the
build if any source isn't listed in the subsys's
`SRCS`.  The top-level `check-srcs-coverage` is a
prerequisite of `build` and dispatches each subsys's
check.  `test/Makefile` carries the analogous
`check-tests-coverage` drift guard for the TESTS
manifest.

#### Adding a new subsystem

1. Create `<subsys>/Makefile` modelled on
   `lib/Makefile` (the simplest existing one).  Set
   `SUBSYS := <subsys>`, list `SRCS`, derive
   `OBJS`, define `all` / `shared` /
   `check-srcs-coverage` / `clean`.
2. In the top-level Makefile's Section 2, add
   `<SUBSYS>OBJS` to the union block.
3. In Section 3, add a phony dispatch target
   `<subsys>:` that invokes `$(MAKE) -C <subsys>`,
   plus the appropriate dep edges to upstream
   subsystems.
4. If the new subsys links a binary, add its build
   rule to that subsys's Makefile (model on
   `ph2run/Makefile`).

`make clean` removes `$(BUILDDIR)/`.  `make
distclean` additionally removes the built binaries
and `libphase2.so`.  Sub-makes refuse to run
standalone -- they error with a pointer to the
top-level Makefile.
