phase2
======

[![CI](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml)

Hamiltonian simulation with Trotter product formula. 🏭


Dependencies
------------

To run a simulation, you will need:

- CPU supporting AVX2 instuction set extension
- Linux x86-64 platform, with a C11/C++11 compiler toolchain
- [OpenMPI][openmpi-website]
- Parallel [HDF5][hdf5-website]

Assuming we're on Ubuntu 22.04 or later, you can install the dependencies with:

```bash
sudo apt install gcc
sudo apt install libopenmpi-dev openmpi-common
sudo apt install libhdf5-dev hdf5-tools libhdf5-mpi-dev libhdf5-openmpi-dev
```

If you want to run the simulation on Euler cluster (ETHZ), you only need to
load necessary modules:

```bash
ml load stack/2024-05
ml load curl/8.4.0-mxgyalo
ml load libszip/2.1.1-xxz4pi5
ml load openmpi
ml load hdf5
ml load python
```

*Note.* This list assumes that Euler has fully migrated to Ubuntu by now
(expected to complete by 8 July 2024). If you log in to the legacy CentOS-7
version, the procedure (and how to compile the sources) is a bit more
complicated.


[hdf5-website]: https://www.hdfgroup.org/solutions/hdf5/
[openmpi-website]: https://www.open-mpi.org/


Compiling the source code
-------------------------

Download and compile the source code with:

```bash
git clone https://github.com/Quantum-for-Life/phase2
cd phase2
make
```

To run the test suite:

```bash
make test
```

Consult [Makefile](./Makefile) for how to configure the build system.


How to use it
-------------

Prepare an input file for the simulation according to the
[specification](./simul/simul-h5-specs.md). Then run:

```bash
mpirun -n [NUM_CPUS] ./ph2run/ph2run-trott [SIMUL_FILE] [NUM_STEPS]
```

or

```bash
mpirun -n [NUM_CPUS] ./ph2run/ph2run-qdrift [SIMUL_FILE]
```

where

- `NUM_CPUS` is the `mpirun` option specifying the number of processes to
  run. This must be a power of 2.
- `SIMUL_FILE` is the path to the simulation file in the HDF5 format.
- `NUM_STEPS` is a integer number of Trotter steps the program is going to
  compute.

Optionally, you can specify the level of log messages for the program to report,
by setting `PHASE2_LOG` environment variable to be one of: `trace`,
`debug`, `info`, `warn`, `error`, `fatal`. E.g.,

```bash
PHASE2_LOG=info mpirun -n 8 ./ph2run-trott simul.h5 100
```

will compute 100 Trotter steps for a Hamiltonian specified in the file
`simul.h5` using 8 MPI processes, and write the result back to the same file.

See also the directory: [./simul](./simul) for an example of a simple automated system,
and the repository:
[Quantum-for-Life/simul-trott-error](https://github.com/Quantum-for-Life/simul-trott-error)
for a real-life application.

