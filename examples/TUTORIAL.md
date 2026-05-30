# Tutorial: ground-state energy with phase2

This walks the full pipeline end to end: turn a molecular integral
file into a simulator input, evolve it under a Trotter circuit, and
read the ground-state energy off the result. It uses the bundled water
CAS(5,6) fixture; the exact ground-state energy is **-74.997 Ha**, so
you know what to expect.

## 0. Build

From the repository root:

```sh
make                            # build/ph2run/ph2run
python3 -m venv ~/.venv && . ~/.venv/bin/activate
pip install -e ".[examples]"    # h5py, numpy, scipy
```

## 1. The input file

A run reads one HDF5 file, `simul.h5`, carrying:

- `/pauli_hamil` — the Hamiltonian as `H = sum_k c_k P_k`, stored as
  normalised coefficients plus two scalars: `normalization`
  (`1/sum|c_k|`, so the simulator evolves a unit-norm operator) and
  `offset` (the identity term + nuclear repulsion).
- `/state_prep/multidet` — the reference state `|psi>`.

The committed `data/hamil.h5` already has both, prepared from
`data/FCIDUMP` (the integrals) and `data/INPUTST` (the reference
determinant). To rebuild it yourself you need the `[prep]` extras
(`qiskit-nature`, `pyscf`):

```sh
make -C examples/simul/trott regen
```

`parse_fcidump.py` applies the Jordan-Wigner transform; the resulting
`/pauli_hamil` is a 10-qubit, 251-term operator.

## 2. Evolve

The simulator computes the overlap of the time-evolved reference state
with itself:

    values[s] = <psi| exp(i * delta * H_norm)^(s+1) |psi>

one complex number per Trotter step, at uniform spacing `dt = delta`.
Run 1st-order Trotter for 4096 steps on two MPI ranks:

```sh
make -C examples/simul/trott
# = cp data/hamil.h5 simul.h5
#   mpirun -n 2 ph2run -S simul.h5 trott --delta=0.1 --steps=4096
```

Results land in `/circ_trott/values` inside `simul.h5`.

## 3. Read the energy: FFT

`values[s]` is a sum of tones, one per eigenstate in `|psi>`:

    values[s] = sum_n |<psi|n>|^2 exp(i * lambda_n * delta * s)

A discrete Fourier transform resolves them. A peak at FFT frequency
`f` corresponds to an eigenvalue `lambda = 2*pi*f / delta` of the
normalised Hamiltonian, hence a physical energy

    E = 2*pi*f / (delta * normalization) + offset

The dominant peak (largest magnitude, ignoring the DC bin) is the
eigenstate with the largest overlap with `|psi>`. For a reference
determinant close to the ground state that peak **is** the
ground-state energy. `energy_fft.py` does exactly this:

```sh
make -C examples/simul/trott analyze     # already run by step 2
# -75.00018314380438
```

`-75.000` against the exact `-74.997` — the residual is the FFT bin
width, `2*pi / (steps * delta * normalization)`. Raise `--steps`
(or `delta`) to land the line closer to a bin centre. See the full
spectrum with:

```sh
python3 examples/scripts/energy_fft.py --peaks examples/simul/trott/simul.h5
```

The 2nd-order method (`make -C examples/simul/trott2`) uses the same
readout on the same time grid and is more accurate per step.

## 4. The randomised algorithms: averaged overlap

qDRIFT and the composite integrator do not produce a time series: each
recorded value is an independent Monte-Carlo *sample* of the overlap
at one fixed effective evolution time `T`. Averaging the samples
recovers the overlap, and its phase gives the energy:

    E = arg(mean(values)) / (T * normalization) + offset

with `T = depth * asin(step_size)` for qDRIFT and `T = steps *
angle_det` for the composite. `energy_mc.py` reads this:

```sh
make -C examples/simul/qdrift     # -2.757,-74.976
make -C examples/simul/cmpsit     # -2.823,-75.043
```

(CSV: normalised eigenvalue, then physical energy.) The estimate is
reliable while the averaged overlap stays coherent — keep `T` of order
1 and use enough samples; statistical error falls as `1/sqrt(samples)`.

## Where things live

| Piece            | Path                              |
|------------------|-----------------------------------|
| simulator        | `build/ph2run/ph2run`             |
| input prep       | `examples/scripts/parse_*.py`     |
| FFT energy       | `examples/scripts/energy_fft.py`  |
| Monte-Carlo energy | `examples/scripts/energy_mc.py` |
| pipelines        | `examples/simul/<algo>/Makefile`  |
| fixture          | `examples/data/`                  |

For the pure-library (no HDF5) interface, see `pauli_rotation.py` and
[../doc/python.md](../doc/python.md).
