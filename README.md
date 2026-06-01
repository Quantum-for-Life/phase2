phase2
======

[![CI](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml)

Full state-vector simulation of Hamiltonian time evolution,
targeting quantum phase estimation workflows in computational
chemistry.  Runs on CPU and NVIDIA GPU clusters via MPI and
HDF5; reference algorithms include 1st- and 2nd-order Trotter,
qDRIFT, and a composite (partially randomised) integrator.

Reference: arXiv:[2504.17881](https://arxiv.org/abs/2504.17881).


Install
-------

Build dependencies on Ubuntu 22.04 or later:

```bash
sudo apt install gcc libopenmpi-dev openmpi-common \
                 libhdf5-dev hdf5-tools
```

On the ETH Euler cluster:

```bash
ml load stack/2024-06 openmpi/4.1.6 hdf5/1.14.3 \
        python/3.11.6 curl/8.4.0-s6dtj75 \
        libszip/2.1.1-gz5ijo3 zlib/1.3-mktm5vz
```

Build:

```bash
git clone https://github.com/Quantum-for-Life/phase2
cd phase2
make
```


Run
---

Simulate four 1st-order Trotter steps on the water
CAS(5,6) test fixture, on two MPI ranks:

```bash
mpirun -n 2 ./build/ph2run/ph2run -S test/data/H2O_CAS56.h5 \
       trott -D 0.1 -s 4
```

Results land in the `/circ_trott` group of the same HDF5
file.  Other subcommands — `trott2`, `qdrift`, `cmpsit` —
share the same `-S FILE` convention; see
`./build/ph2run/ph2run --help` and `./build/ph2run/ph2run CMD --help`
for the flag surface.

For scripted prepare-run-analyse pipelines (FCIDUMP to energy
estimate) for each algorithm, see [examples/simul](examples/simul).

The MPI rank count must be a power of two and must not
exceed `2^(nqb-1)`, where `nqb` is the qubit width of the
register.  Default log verbosity is `info` on stdout;
raise with `PHASE2_LOG=debug` (requires `make debug` for
trace/debug to compile in), enable per-rank logging in
distributed runs with `PHASE2_LOG_ALL=1`.  See
[doc/logging.md](doc/logging.md) for the full reference.


State preparation contract
--------------------------

`simul.h5` carries exactly one initial-state encoding under
`/state_prep`:

- `/state_prep/multidet` — explicit list of (bitstring,
  complex amplitude) pairs.
- `/state_prep/coeff_matrix` — real `(n_sites, n_occ)`
  molecular-orbital coefficient matrix.  The simulator
  expands it on the fly via Slater-Condon, scattering
  amplitudes directly into the MPI-distributed register.
  Supports closed/open shell, Z₂ + S_z tapering, and CSF
  superpositions.

The file is rejected at open time if both subgroups are
present, or if neither is.  See
[doc/simul-h5-specs.md](doc/simul-h5-specs.md) for the full
schema, dispatch rules, and tapering convention.


GPU build
---------

Build against CUDA-aware OpenMPI:

```bash
make clean
make BACKEND=cuda
make BACKEND=cuda check
```

Requires one GPU per MPI process and a CUDA-aware MPI
installation.  Override `CUDA_PREFIX` if the toolkit lives
outside `/usr/local/cuda`.


Python
------

`phase2` ships a thin Python wrapper around a single
overlap-evaluation entry point:

```bash
make shared
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
```

```python
import phase2
r = phase2.run(["X0", "Z1"], [0.3, -0.2], 1.0, "00")
```

With MPI (ranks must be a power of two and at most half the
number of amplitudes):

```bash
mpirun -n 4 python examples/pauli_rotation.py
```

See [doc/python.md](doc/python.md) for the full API
reference and limitations.


Tests and benchmarks
--------------------

```bash
make check                 # full test suite (parallel, cargo-style)
make check-mpi MPIRANKS=4  # same, each test under mpirun -n 4
make check-data            # subset filter (e.g. check-paulis, check-circ)
make bench                 # paulis and qreg micro-benchmarks
make test-asan             # ASan + UBSan build
make test-valgrind         # leak/UB check via valgrind
```

Test fixtures live in `test/data/`.  See
[doc/testing.md](doc/testing.md) for the harness
contract, runner flags, and recipes for adding tests.


Documentation
-------------

- [doc/phase2.md](doc/phase2.md) — reference manual:
  design, kernels (CPU and CUDA), algorithms, API, build.
- [doc/logging.md](doc/logging.md) — logging subsystem:
  format, levels, env vars, MPI, policy for new code.
- [doc/testing.md](doc/testing.md) — test subsystem:
  harness, parallel runner, conventions, recipes.
- [doc/python.md](doc/python.md) — Python interface.
- [doc/simul-h5-specs.md](doc/simul-h5-specs.md) — HDF5
  input/output schema.


Citation
--------

If you use phase2 in published work, cite:

> Marek Miller, Jakob Gunther, Freek Witteveen, Matthew S.
> Teynor, Mihael Erakovic, Markus Reiher, Gemma C. Solomon,
> and Matthias Christandl, *phase2: Full-State Vector
> Simulation of Quantum Time Evolution at Scale*,
> arXiv:2504.17881.


Maintainer & licence
--------------------

Marek Miller <mlm@math.ku.dk>.  BSD 3-Clause — see
[LICENSE](./LICENSE).  Report bugs and submit patches via
[GitHub Issues](https://github.com/Quantum-for-Life/phase2/issues).
