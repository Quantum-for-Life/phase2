# phase2 examples

Two kinds of example live here:

- **`pauli_rotation.py`** — a self-contained demo of the
  `libphase2.so` Python (ctypes) interface.  It evaluates
  overlaps of Pauli-rotation products against computational
  basis states and checks each result against an analytic
  value.  No HDF5, no input files.

- **`simul/`** — end-to-end simulation pipelines driven by
  `ph2run`.  Each algorithm (`trott`, `trott2`, `qdrift`,
  `cmpsit`) has a `Makefile` that prepares `simul.h5`, runs
  the simulation under MPI, and estimates the ground-state
  energy with a Python analysis script.  See
  [simul/README.md](simul/README.md).

Shared assets:

- `data/` — the molecular fixture: `FCIDUMP` + `INPUTST`
  (source of truth) and the prepared `hamil.h5`
  (`/pauli_hamil` + `/state_prep/multidet`).
- `scripts/` — input preparation (`parse_fcidump.py`,
  `parse_inputst.py`) and analysis (`trott_rpe.py`,
  `trott_fft.py`, `qdrift_rpe.py`, `qdrift_sample.py`,
  `cmpsit_rpe.py`).

## Prerequisites

Build the simulator and (for the ctypes demo) the shared
library from the repository root:

```sh
make             # build/ph2run/ph2run
make shared      # build/libphase2.so   (pauli_rotation.py only)
```

Python dependencies (use the project virtualenv):

```sh
pip install -e ".[examples]"   # h5py, numpy, scipy: run + analyse
pip install -e ".[prep]"       # + qiskit-nature, pyscf: regenerate hamil.h5
```

## Quick start

```sh
python examples/pauli_rotation.py          # library demo
make -C examples/simul/trott               # Trotter pipeline
```

The `simul/` pipelines need only the `[examples]` extras and
the committed `hamil.h5`; `qiskit-nature`/`pyscf` are required
only to rebuild `hamil.h5` from `FCIDUMP` (`make regen`).
