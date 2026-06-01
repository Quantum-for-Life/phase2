# `simul` pipelines

End-to-end ground-state energy estimation with `ph2run`. See
[../TUTORIAL.md](../TUTORIAL.md) for the narrated walkthrough. Each
subdirectory drives one algorithm through the same three stages:

1. **Prepare** — seed `simul.h5` from the shared
   `../data/hamil.h5` (`/pauli_hamil` + `/state_prep`).
2. **Simulate** — run `ph2run <algo>` under `mpirun`.
3. **Analyse** — extract the energy from `/circ_<algo>/values`.

Energy extraction uses one method per data shape:

| Directory | Algorithm                     | Analysis        |
|-----------|-------------------------------|-----------------|
| `trott`   | 1st-order Trotter             | `energy_fft.py` |
| `trott2`  | 2nd-order (Strang) Trotter    | `energy_fft.py` |
| `qdrift`  | qDRIFT (randomised)           | `energy_mc.py`  |
| `cmpsit`  | composite (det. + randomised) | `energy_mc.py`  |

`trott`/`trott2` record a uniform overlap time series, so the energy
is the dominant peak of its FFT. `qdrift`/`cmpsit` record independent
Monte-Carlo samples at one effective time, so the energy is the phase
of their averaged overlap.

## Running

From the repository root, build `ph2run` and install the `[examples]`
extras (see [../README.md](../README.md)), then:

```sh
make -C examples/simul/trott                     # default parameters
make -C examples/simul/trott  DELTA=0.1 TROTT_STEPS=4096 MPI_RANKS=2
make -C examples/simul/qdrift STEP_SIZE=0.0625 DEPTH=16 SAMPLES=1024
make -C examples/simul/cmpsit LENGTH=200 STEPS=8 SAMPLES=256
make -C examples/simul/<algo> clean
```

Each pipeline writes its energy to `<algo>/simul.h5.proc` (also echoed
to the terminal). The MPI rank count must be a power of two and at
most `2^(nqb-1)`; the bundled fixture has `nqb = 10`, so `MPI_RANKS`
up to 256 is valid.

## Energy convention

The simulator evolves the **normalised, identity-removed**
Hamiltonian; the physical molecular energy is the recovered eigenvalue
plus `offset`, the scalar term recorded by `parse_fcidump.py`. For the
bundled water CAS(5,6) fixture the exact values are: ground state
**-74.997 Ha**, reference state **-74.963 Ha**.

Both scripts print one CSV line, `E,E_ref,dE`: the estimated energy,
the trial-state energy `E_ref = <psi|H|psi>` (computed directly from
the input by `energy_ref.py`), and `dE = E - E_ref`, the shift the
evolution produced relative to the reference state.

- `energy_fft.py` reports the dominant-peak energy as `E`. Resolution
  is one FFT bin = `2*pi / (steps * delta * norm)`; raise `TROTT_STEPS`
  to sharpen it. `--peaks` lists every detected line.
- `energy_mc.py` reports the averaged-overlap energy as `E`. It is
  reliable while the effective time `T` keeps the mean coherent
  (|mean| not too small) and the per-sample phase stays in `(-pi, pi)`;
  large `T` (or, for `cmpsit`, `angle_rand` far from `angle_det`)
  decoheres the average and biases the readout. Raise `SAMPLES` to
  shrink the `1/sqrt(N)` statistical error.

## Input fixtures

All four pipelines share one molecule, `../data/FCIDUMP` (water,
CAS(5,6): `NORB=5`, `NELEC=6`) and the reference state
`../data/INPUTST`. `make regen` rebuilds `../data/hamil.h5` from them
(needs the `[prep]` extras).

### `FCIDUMP`

Standard FCIDUMP one- and two-electron integrals. `parse_fcidump.py`
maps it to a Pauli Hamiltonian via the Jordan-Wigner transform
(qiskit-nature), storing the normalisation and the scalar `offset`
separately.

### `INPUTST` (input state)

Plain text; each line is one Slater determinant of a convex linear
combination:

```text
F F I I I ... I
```

`F F` are the real and imaginary parts of the complex coefficient; the
`I` values are spin-orbital occupations (`0`/`1`). Alpha and beta
orbitals interleave:

```text
I_{1,alpha} I_{1,beta} I_{2,alpha} I_{2,beta} ...
```

The occupation count must match across rows; whitespace and blank
lines are ignored.
