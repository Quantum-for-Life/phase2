# phase2

[![CI](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml)

Simulate generic Hamiltonians using Trotter product formula. 🏭

# Dependencies

To run the simulation, you will need:

- Linux platform, with a C/C++ compiler toolchain
- [OpenMPI][openmpi-website]
- Parallel [HDF5][hdf5-website]

Assuming we're on x86-64 Ubuntu, you can install the dependencies with:

```bash
sudo apt install gcc
sudo apt install libopenmpi-dev openmpi-common
sudo apt install libhdf5-dev hdf5-tools libhdf5-mpi-dev libhdf5-openmpi-dev 
```

If you want to run the simulation on Euler cluster (ETHZ), you only need to
load the necessary modules:

```bash
env2lmod
module load gcc
module load hdf5
module load openmpi
```

[hdf5-website]: https://www.hdfgroup.org/solutions/hdf5/

[openmpi-website]: https://www.open-mpi.org/

# Compiling the source code

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

Consult [Makefile](./Makefile) for how to configure the build system, if you
have the dependencies install in a different location.

# How to use it

Prepare an input file for the simulation according to the
[specification](./simul/simul-h5-specs.md). Then run:

```bash
mpirun -n [NUM_CPUS] ./ph2run/ph2run-trott [SIMUL_FILE] [NUM_STEPS]
```

where

- `NUM_CPUS` is the `mpirun` option specifying the number of processes to
  run. This must be a power of 2.
- `SIMUL_FILE` is the path to the simulation file in the HDF5 format
- `NUM_STEPS` is a integer number of Trotter steps the program is going to
  compute.

Optionally, you can specify the level of log messages for the program to report,
by setting `PHASE2_LOG` environment variable to be one of: `trace`,
`debug`, `info`, `warn`, `error`, `fatal`. E.g.,

```bash
PHASE2_LOG=info mpirun -n 8 ./ph2run simul.h5 100 
```

will compute 100 Trotter steps for a Hamiltonian specified in the file
`simul.h5` using 8 MPI processes, and write the result to the same file.

See also the directory: `./simul` for an example of a simple automated system.

*TODO*:

- Memory requirements
- Single precision floating point calculation is a CMake switch: `-DPRECISION=1`
