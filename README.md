# phase2

[![CI](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml)

Quantum phase estimation (QPE) with parallel computing. 🏭

# Dependencies

To run the simulation, you will need:

- Linux platform, with a C/C++ compiler toolchain
- CMake>=3.9
- [HDF5][hdf5-website]
- optional: [OpenMPI][openmpi-website]

To prepare input for the simulation and to run the test suite:

- Python >= 3.10 and a virtual environment with several packages installed
  locally (see [README.md](./python/README.md) in the `python/` subsystem for
  more details.

Assuming we're on Ubuntu, you can install the dependencies in one go:

```bash
sudo apt install gcc g++ cmake 
sudo apt install libhdf5-dev hdf5-tools
```

and optionally:

```bash
sudo apt install libopenmpi-dev openmpi-common
sudo apt install libhdf5-mpi-dev libhdf5-openmpi-dev 
```

If you want to run the simulation on Euler cluster (ETHZ), you only need to load
necessary modules:

```bash
env2lmod
module load gcc
module load cmake
module load hdf5
module load openmpi
```

[hdf5-website]: https://www.hdfgroup.org/solutions/hdf5/

[openmpi-website]: https://www.open-mpi.org/

## Getting the sources

You can clone the repo and its dependencies:

```bash
git clone --recurse-submodules https://github.com/Quantum-for-Life/phase2.git
```

If you prefer to access GitHub via SSH instead, be aware the repositories
this package depends on as submodules have their URLs hard-coded as HTTPS.
And some of them are still private to _Quantum-for-Life_.

# Compiling the source code

```bash
cd phase2
mkdir build && cd build
cmake ..
make
```

If you want to enable MPI support, set `DISTRIBUTED` flag:

```bash
cmake -DDISTRIBUTED=ON ..
make
```

# Testing

Run the test suite by typing:

```bash
cmake -DBUILD_TESTING=ON ..
make && ctest -V
```

# How to use it

[TBA]
