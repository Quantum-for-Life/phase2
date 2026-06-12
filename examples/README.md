# phase2 examples

- **[TUTORIAL.md](TUTORIAL.md)** — start here: FCIDUMP → simulate →
  ground-state energy, end to end on the bundled fixture.

- **`pauli_rotation.py`** — self-contained demo of the `libphase2.so`
  Python (ctypes) interface: overlaps of Pauli-rotation products
  against basis states, each checked analytically. No HDF5, no input
  files. See [../doc/python.md](../doc/python.md).

- **`simul/`** — `ph2run` simulation pipelines, one `Makefile` per
  algorithm (`trott`, `trott2`, `qdrift`, `cmpsit`). See
  [simul/README.md](simul/README.md).

Shared assets:

- `data/` — the molecular fixture: `FCIDUMP` + `INPUTST` (source of
  truth) and the prepared `hamil.h5` (`/pauli_hamil` +
  `/state_prep/multidet`).
- input preparation and energy analysis go through the worksheet
  toolkit `util/ph2.py` (`ph2 hamil`, `ph2 stprep`,
  `ph2 energy {fft,mc,ref,rpe}`); see [../doc/ph2.md](../doc/ph2.md).

## Prerequisites

```sh
make             # build/ph2run/ph2run
make shared      # build/libphase2.so   (pauli_rotation.py only)

pip install -e ".[examples]"   # h5py, numpy, scipy: run + analyse
pip install -e ".[prep]"       # + qiskit-nature, pyscf: regenerate hamil.h5
```

## Quick start

```sh
python examples/pauli_rotation.py     # library demo
make -C examples/simul/trott          # Trotter pipeline -> energy
```

The `simul/` pipelines need only the `[examples]` extras and the
committed `hamil.h5`; `qiskit-nature`/`pyscf` are required only to
rebuild `hamil.h5` from `FCIDUMP` (`make regen`).
