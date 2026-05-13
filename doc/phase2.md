# phase2

Full-state vector quantum Hamiltonian simulation library for
quantum phase estimation.

Version 0.12.1

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
                     libhdf5-openmpi-dev

The parallel HDF5 build (`libhdf5-openmpi-dev`) is
essential; the serial variant will not work.

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

### 3.5 Backend Abstraction

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

**`PHASE2_LOG`**: Controls the logging verbosity.  Accepted
values (case-sensitive): `trace`, `debug`, `info`, `warn`,
`error`, `fatal`.  The log level is read once during
`world_init`.  Default behaviour: info-level messages are
emitted.

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

### 8.4 Data I/O (`include/phase2/data.h`)

```c
typedef int64_t data_id;
#define DATA_INVALID_FID (-1)
```

---

```c
data_id data_open(const char *filename);
```

Open an HDF5 simulation file for reading and writing.
Collective MPI operation (all ranks must call).

**Return value:** A valid file handle, or DATA_INVALID_FID
on error.

---

```c
void data_close(data_id fid);
```

Close a previously opened HDF5 file.

---

```c
int data_grp_create(data_id fid, const char *grp_name);
```

Create an HDF5 group in the root of the file.

**Return value:** 0 on success, -1 on error.

---

```c
int data_attr_read_i(data_id fid, const char *grp_name,
    const char *attr_name, int *a);
int data_attr_read_ul(data_id fid, const char *grp_name,
    const char *attr_name, unsigned long *a);
int data_attr_read_dbl(data_id fid, const char *grp_name,
    const char *attr_name, double *a);
```

Read a scalar attribute from the named group.  A type-
generic macro `data_attr_read` dispatches based on the
pointer type of the fourth argument.

**Return value:** 0 on success, -1 on error.

---

```c
int data_attr_write_i(data_id fid, const char *grp_name,
    const char *attr_name, int a);
int data_attr_write_ul(data_id fid, const char *grp_name,
    const char *attr_name, unsigned long a);
int data_attr_write_dbl(data_id fid, const char *grp_name,
    const char *attr_name, double a);
```

Create and write a scalar attribute in the named group.
A type-generic macro `data_attr_write` dispatches based on
the type of the fourth argument.

**Return value:** 0 on success, -1 on error.

---

```c
int data_multidet_getnums(data_id fid, uint32_t *nqb,
                          size_t *ndets);
```

Retrieve the number of qubits and determinants from the
`/state_prep/multidet` group.

**Return value:** 0 on success, -1 on error.

---

```c
int data_multidet_foreach(data_id fid,
    int (*op)(_Complex double cf, uint64_t idx, void *),
    void *op_data);
```

Iterate over each determinant in the multidet group.  For
each determinant, calls `op(cf, idx, op_data)` where `cf`
is the complex coefficient and `idx` is the basis state
index.

**Return value:** 0 if the full iteration completed, -1 on
data retrieval error, or the nonzero return value of `op`
if iteration was terminated early.

---

```c
int data_hamil_getnums(data_id fid, uint32_t *nqb,
                       size_t *nterms);
```

Retrieve the number of qubits and Hamiltonian terms from
the `/pauli_hamil` group.

**Return value:** 0 on success, -1 on error.

---

```c
int data_hamil_getnorm(data_id fid, double *norm);
```

Retrieve the normalisation factor from the `/pauli_hamil`
group.  All Hamiltonian coefficients are multiplied by this
factor during loading.

**Return value:** 0 on success, -1 on error.

---

```c
int data_hamil_foreach(data_id fid,
    int (*op)(double, unsigned char *, void *),
    void *op_data);
```

Iterate over each Hamiltonian term.  For each term, calls
`op(coeff, paulis_array, op_data)` where `coeff` is the
real coefficient and `paulis_array` is a temporary array of
`num_qubits` bytes encoding single-qubit Paulis (0=I, 1=X,
2=Y, 3=Z).

**Return value:** 0 if the full iteration completed, -1 on
error, or the nonzero return value of `op` if terminated
early.  The `paulis_array` pointer is invalidated after
iteration.

---

```c
int data_res_write(data_id fid, const char *grp_name,
    const char *dset_name, const _Complex double *vals,
    size_t nvals);
```

Write an array of complex values as a dataset in the named
group.  The dataset has shape (nvals, 2) with columns for
real and imaginary parts.

**Return value:** 0 on success, -1 on error.

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
};
```

Complete simulation context: Hamiltonian, initial state,
quantum register, cache, and output values buffer.

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
int circ_init(struct circ *ct, data_id fid,
              size_t vals_len);
```

Initialise a complete simulation context from an HDF5 file.
Loads the Hamiltonian and multideterminant state, allocates
the quantum register, initialises the cache, and allocates
the output buffer.

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
struct trott { struct circ ct; struct trott_data dt; };
```

---

```c
int trott_init(struct trott *tt,
    const struct trott_data *dt, data_id fid);
```

Initialise a Trotter simulation from HDF5 data.  Loads
Hamiltonian and initial state, sorts the Hamiltonian
lexicographically.

**Return value:** 0 on success, -1 on error.

---

```c
void trott_free(struct trott *tt);
```

Free all Trotter simulation resources.

---

```c
int trott_simul(struct trott *tt);
```

Run the Trotter simulation for `dt.steps` steps,
populating `ct.vals` with overlap measurements.

**Return value:** 0 on success, -1 on error.

---

```c
int trott_write_res(struct trott *tt, data_id fid);
```

Write results to the `/circ_trott` group in the HDF5 file.

**Return value:** 0 on success, -1 on error.

---

### 8.8 Symmetric Trotter (`include/circ/trott2.h`)

Strang 2nd-order product formula.  One step is a forward
sweep at `delta/2` followed by a reverse sweep at
`delta/2`, both expressed through `circ_step` /
`circ_step_reverse`.

```c
struct trott2_data { double delta; size_t steps; };
struct trott2 { struct circ ct; struct trott2_data dt; };

int  trott2_init(struct trott2 *t2,
        const struct trott2_data *dt, data_id fid);
void trott2_free(struct trott2 *t2);
int  trott2_simul(struct trott2 *t2);
int  trott2_write_res(struct trott2 *t2, data_id fid);
```

Behaviour mirrors `trott_*` with results written to
`/circ_trott2`.  See `include/circ/trott2.h` for the
per-function contract.

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
};
```

---

```c
int qdrift_init(struct qdrift *qd,
    const struct qdrift_data *dt, data_id fid);
```

Initialise a qDRIFT simulation.  Loads Hamiltonian, builds
the CDF, initialises the PRNG.

**Return value:** 0 on success, -1 on error.

---

```c
void qdrift_free(struct qdrift *qd);
```

Free all qDRIFT simulation resources.

---

```c
int qdrift_simul(struct qdrift *qd);
```

Run the qDRIFT simulation for `dt.samples` independent
samples, each of depth `dt.depth`.

**Return value:** 0 on success, -1 on error.

---

```c
int qdrift_write_res(struct qdrift *qd, data_id fid);
```

Write results to the `/circ_qdrift` group in the HDF5
file.

**Return value:** 0 on success, -1 on error.

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
};
```

---

```c
int cmpsit_init(struct cmpsit *cp,
    const struct cmpsit_data *dt, data_id fid);
```

Initialise a composite simulation.  Loads Hamiltonian,
splits into deterministic and randomised parts, builds the
CDF for the randomised part, sorts the deterministic part
lexicographically, initialises the PRNG.

**Return value:** 0 on success, -1 on error.

---

```c
void cmpsit_free(struct cmpsit *cp);
```

Free all composite simulation resources.

---

```c
int cmpsit_simul(struct cmpsit *cp);
```

Run the composite simulation: for each of `dt.samples`
samples, apply `dt.steps` second-order Trotter steps and
measure.

**Return value:** 0 on success, -1 on error.

---

```c
int cmpsit_write_res(struct cmpsit *cp, data_id fid);
```

Write results to the `/circ_cmpsit` group in the HDF5
file.

**Return value:** 0 on success, -1 on error.

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
happens during data loading in `circ_hamil_from_file`.

---

## 10. Building and Testing

### 10.1 Make Targets

| Target        | Description                           |
|---------------|---------------------------------------|
| `all`         | Build everything (programs, benchmarks, tests) |
| `build`       | Build `ph2run` and dependencies       |
| `build-test`  | Build the test suite                  |
| `build-bench` | Build benchmarks                      |
| `check`       | Run all tests (single rank)           |
| `check-mpi`   | Run all tests under MPI               |
| `bench`       | Run benchmarks (single rank)          |
| `bench-mpi`   | Run benchmarks under MPI              |
| `debug`       | Build with debug flags (-g -Og)       |
| `clean`       | Remove object files                   |
| `distclean`   | Remove object files, binaries, tests  |
| `format`      | Run clang-format on all sources       |

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

The Python harness `test/t-ref-coeff_matrix.py` cross-checks
the C expansion against an in-tree reference oracle in
`test/ref/coeff_matrix_reference.py`; it is run by
`make check-python`.

Test data files reside in `test/data/`.  The slow test
`t-state_prep_coeff_large` is built and run separately via
`make build-test-slow` / `make check-slow`.  Sanitiser and
valgrind variants are available as `make test-asan`,
`make test-valgrind`, `make test-mpi-asan`.  Tests are
built with `-DDEBUG -g -Og` flags.
