# phase2 Python Interface

## Overview

The `phase2` Python module provides access to the phase2
Hamiltonian simulation engine.  It wraps the C function
`phase2_run` from `libphase2.so` via `ctypes`, exposing a
single function `phase2.run()` that computes overlaps of
Trotterised Pauli rotation products with computational
basis states.

The core computation is:

```
<psi| prod_k exp(i * delta * coeffs[k] * P_k) |psi>
```

where `P_k` are Pauli strings, `coeffs[k]` are real
coefficients, `delta` is a step-size parameter, and `|psi>`
is a computational basis state.  The product is ordered
left-to-right: `k=0` is applied first (innermost), `k=last`
is applied last (outermost).

All compute-intensive work runs in C with optional MPI
parallelism.  Python serves as the high-level orchestration
layer.

## Installation

### 1. Build the shared library

```sh
make shared
```

This produces `libphase2.so` under `build/`.  The build
requires GCC, OpenMPI, and parallel HDF5 development headers.

### 2. Create and activate a virtual environment

```sh
python3 -m venv .venv
source .venv/bin/activate
```

The virtual environment is local to the project directory
(`.venv/`).  Activate it in each new shell session with
`source .venv/bin/activate`.

### 3. Install the Python package

```sh
pip install .
```

Or, for development (editable install):

```sh
pip install -e .
```

To install with test dependencies:

```sh
pip install -e ".[test]"
```

### Verify the installation

```sh
pytest -v
```

The test suite is at `python/tests/test_phase2.py`.  The
`testpaths` setting in `pyproject.toml` ensures `pytest`
finds it automatically.

### Complete setup sequence

```sh
make shared
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
pytest -v
```

### 4. Library discovery

The module locates `libphase2.so` by searching, in order:

1. The path in the `PHASE2_LIB` environment variable.
2. `../../build/libphase2.so` relative to the installed
   package -- the in-tree development build's output.
3. `../../libphase2.so` relative to the installed package.
4. `../libphase2.so` relative to the installed package.
5. `./libphase2.so` relative to the installed package.

To override the search, set `PHASE2_LIB` explicitly:

```sh
export PHASE2_LIB=/path/to/libphase2.so
```

## Quick Start

```python
import phase2

result = phase2.run(
    paulis=["Z0 Z1", "X0", "X1"],
    coeffs=[0.5, -0.3, -0.3],
    delta=0.01,
    psi="00",
)
print(result)  # complex number
```

## API Reference

### `phase2.run(paulis, coeffs, delta, psi) -> complex`

Compute the overlap of a product of Pauli rotations with a
computational basis state.

**Parameters:**

- `paulis` (`list[str]`): List of Pauli strings.  Each
  string specifies a Pauli operator as space-separated
  tokens (see Pauli String Format below).  The number of
  terms is `len(paulis)`.

- `coeffs` (`list[float]`): Real coefficients, one per
  Pauli string.  Must satisfy `len(coeffs) == len(paulis)`.

- `delta` (`float`): Scaling factor (time step).  The
  rotation angle for term `k` is `delta * coeffs[k]`.

- `psi` (`str`): Computational basis state as a bit string
  of `'0'` and `'1'` characters.  The number of qubits is
  `len(psi)`.  For example, `"000"` is the 3-qubit ground
  state `|000>`.

**Returns:**

- `complex`: The overlap value
  `<psi| prod_k exp(i * delta * coeffs[k] * P_k) |psi>`.

**Raises:**

- `ValueError`: If `len(paulis) != len(coeffs)`.
- `ValueError`: If the C library returns an error code
  (invalid input or allocation failure).

**Product ordering:**

The product applies terms left to right.  Term `k=0` acts
on `|psi>` first (innermost), and the last term acts last
(outermost).  For non-commuting operators, this ordering
matters.

## Pauli String Format

Each Pauli string is a space-separated sequence of tokens.
Each token is a Pauli operator (`X`, `Y`, or `Z`) followed
immediately by a zero-indexed qubit number.

Qubits not mentioned in the string are treated as identity.

**Examples:**

| String       | Operator                           |
|--------------|------------------------------------|
| `"Z0"`       | Z on qubit 0                       |
| `"X0"`       | X on qubit 0                       |
| `"Z0 Z1"`   | Z tensor Z on qubits 0 and 1      |
| `"X0 Z1 Y2"`| X tensor Z tensor Y on qubits 0-2 |
| `"X5"`       | X on qubit 5, identity elsewhere   |

**Rules:**

- Qubit indices are non-negative integers.
- The number of qubits is determined by `len(psi)`, not
  by the indices in the Pauli strings.
- Each qubit may appear at most once per string.
- Valid operators are `X`, `Y`, and `Z`.  The identity `I`
  is implicit for unmentioned qubits.

## MPI Usage

The C library supports MPI parallelism.  When invoked via
`mpirun`, the state vector is distributed across ranks.
MPI is initialised automatically on the first call to
`phase2.run()`.

### Single-rank execution

```sh
python examples/pauli_rotation.py
```

No special configuration is required.  The library detects
that MPI is not active and runs on a single process.

### Multi-rank execution

```sh
mpirun -n 2 python examples/pauli_rotation.py
mpirun -n 4 python examples/pauli_rotation.py
```

**Constraints:**

- The number of MPI ranks must be a power of two.
- The number of ranks must not exceed `2^(nqb-1)`, where
  `nqb` is the number of qubits (`len(psi)`).  Each rank
  must hold at least 2 amplitudes (1 lo-qubit).  For
  example, `mpirun -n 2` requires at least 2 qubits.

### Example with environment variable

```sh
export PHASE2_LIB=/opt/phase2/lib/libphase2.so
mpirun -n 8 python my_simulation.py
```

## Limitations

- **Maximum qubits**: The state vector is stored in memory
  as `2^nqb` complex doubles.  Practical limits depend on
  available RAM.  Single-node simulations are typically
  feasible up to around 30 qubits.  The hard upper bound
  is 64 qubits (`QREG_MAX_WIDTH`).

- **Initial state**: The initial state `|psi>` must be a
  computational basis state specified as a bit string.
  Arbitrary superpositions are not supported as input.

- **Observable**: The function computes the overlap
  `<psi|U|psi>` only.  It does not return the full state
  vector or expectation values of arbitrary observables.

- **Real coefficients**: The coefficients array must contain
  real numbers.  Complex Hamiltonian coefficients are not
  supported.

## Examples

The `examples/` directory contains self-contained scripts:

- `examples/pauli_rotation.py` -- Demonstrates single-qubit
  Z and X rotations, two-qubit ZZ rotations, non-commuting
  Trotter steps, and multi-qubit product rotations.  Each
  section derives the expected result analytically, runs
  `phase2.run()`, and verifies the output.

- `examples/simul/` -- scripted `ph2run` pipelines (Trotter,
  qDRIFT, composite) that prepare `simul.h5`, run under MPI,
  and estimate a ground-state energy.  These use the C
  driver and HDF5 I/O rather than the ctypes wrapper; the
  analysis scripts need the `[examples]` extras (`h5py`,
  `numpy`, `scipy`), and regenerating the input fixture needs
  `[prep]` (`qiskit-nature`, `pyscf`).  See
  [../examples/simul/README.md](../examples/simul/README.md).

