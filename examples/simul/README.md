# `simul` pipelines

End-to-end ground-state energy estimation with `ph2run`.
Each subdirectory drives one algorithm through the same
three stages:

1. **Prepare** — seed `simul.h5` from the shared
   `../data/hamil.h5` (`/pauli_hamil` + `/state_prep`).
2. **Simulate** — run `ph2run <algo>` under `mpirun`.
3. **Analyse** — extract an energy estimate from
   `/circ_<algo>/values` with a script in `../scripts`.

| Directory | Algorithm                     | Analysis        |
|-----------|-------------------------------|-----------------|
| `trott`   | 1st-order Trotter             | `trott_rpe.py`  |
| `trott2`  | 2nd-order (Strang) Trotter    | `trott_rpe.py`  |
| `qdrift`  | qDRIFT (RPE depth sweep)      | `qdrift_rpe.py` |
| `cmpsit`  | composite (det. + randomised) | `cmpsit_rpe.py` |

## Running

From the repository root, build `ph2run` and install the
analysis extras (see [../README.md](../README.md)), then:

```sh
make -C examples/simul/trott                 # default parameters
make -C examples/simul/trott  DELTA=0.05 TROTT_STEPS=512 MPI_RANKS=2
make -C examples/simul/qdrift DELTA=1.0 EPSILON=0.125 SAMPLES=64
make -C examples/simul/cmpsit LENGTH=64 STEPS=8 SAMPLES=128
make -C examples/simul/<algo> clean
```

Each pipeline writes its energy estimate to
`<algo>/simul.h5.proc` (also echoed to the terminal).  The
MPI rank count must be a power of two and at most
`2^(nqb-1)`; the bundled fixture has `nqb = 10`, so
`MPI_RANKS` up to 256 is valid.

## Energy convention

The analysis scripts print CSV.  The simulator evolves the
**normalised, identity-removed** Hamiltonian, so the phase
estimate `E0` is that Hamiltonian's eigenvalue; the total
molecular energy is `E0 + offset`, where `offset` is the
scalar term recorded by `parse_fcidump.py`.  For the bundled
water CAS(5,6) fixture the ground state is near
**-74.96 Ha** (exact diagonalisation: reference state
-74.963, ground state -74.997).

- `trott_rpe.py` / `qdrift_rpe.py` print `E0` and
  `E0 + offset`.
- `cmpsit_rpe.py` prints `E0,E0+offset`.  The composite
  estimate is a single-point phase read and is reliable only
  when the randomised tail shares the deterministic time
  scale (`angle_rand` close to `angle_det`) and the total
  phase stays within `(-pi, pi)`; otherwise it is biased by
  the over-rotated randomised part.  Increase `LENGTH`
  (more deterministic terms) for a cleaner estimate.

## Input fixtures

All four pipelines share one molecule, `../data/FCIDUMP`
(water, CAS(5,6): `NORB=5`, `NELEC=6`) and the reference
state `../data/INPUTST`.  `make regen` rebuilds
`../data/hamil.h5` from them (needs the `[prep]` extras).

### `FCIDUMP`

Standard FCIDUMP one- and two-electron integrals.
`parse_fcidump.py` maps it to a Pauli Hamiltonian via the
Jordan-Wigner transform (qiskit-nature), storing the
normalisation and the scalar `offset` separately.

### `INPUTST` (input state)

Plain text; each line is one Slater determinant of a convex
linear combination:

```text
F F I I I ... I
```

`F F` are the real and imaginary parts of the complex
coefficient; the `I` values are spin-orbital occupations
(`0`/`1`).  Alpha and beta orbitals interleave:

```text
I_{1,alpha} I_{1,beta} I_{2,alpha} I_{2,beta} ...
```

The occupation count must match across rows; whitespace and
blank lines are ignored.
