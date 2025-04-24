phase2
======

[![CI](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml)

Hamiltonian simulation on CPU and GPU clusters. 🏭


Dependencies
------------

To run the simulation, you will need:

- CPU supporting AVX2 instruction set extension
- Linux x86-64 platform, with a C11/C++11 compiler toolchain
- [OpenMPI][openmpi-website]
- Parallel [HDF5][hdf5-website]

Assuming we're on Ubuntu 22.04 or later, you can install the dependencies with:

```bash
sudo apt install gcc
sudo apt install libopenmpi-dev openmpi-common
sudo apt install libhdf5-dev hdf5-tools libhdf5-mpi-dev libhdf5-openmpi-dev
```

If you want to run the simulation on Euler cluster (ETHZ), just load these
modules:

```bash
ml load stack/2024-06
ml load openmpi/4.1.6
ml load hdf5/1.14.3
ml load python/3.11.6
ml load curl/8.4.0-s6dtj75
ml load libszip/2.1.1-gz5ijo3
ml load zlib/1.3-mktm5vz
```

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
make check
```

Individual tests can be found in `./test` directory. The compiled applications
are MPI-aware and can be run in a distributed mode with `mpirun` as well:

```bash
make check-mpi
```

You can specify the number of MPI processes with `MPIRANKS=n` like this:

```bash
make check-mpi MPIRANKS=16
```

The default value is `MPIRANKS=2`. This number *must be a power of two*
and must not exceed the number of cores available.

The system can be configured to use the simulator library
[QuEST](https://github.com/QuEST-Kit/QuEST) as a backend, instead of the
internal engine.  Make sure QuEST is compiled in the `DISTRIBUTED` mode.
Consult [Makefile](./Makefile) for how to configure the build system.

You can also use the software to perform distributed GPU simulation.  You
will need a CUDA-aware MPI implementation on your system.  Follow 
[this tutorial](https://github.com/Quantum-for-Life/cuda-aware-mpi-cluster)
to see how to set up a CUDA-aware MPI cluster.

To use GPUs instead of CPUs, recompile the sources:

```bash
make clean
make BACKEND=cuda
make BACKEND=cuda check
```

The software assumes there is one GPU per MPI process available on local nodes.

On Euler, load the `cuda` module along with those specified above:

```bash
ml load cuda/12.1.1
```


How to use it
-------------

Prepare an input file for the simulation according to the
[specification](doc/simul-h5-specs.md). Then run:

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
- `NUM_STEPS` is an integer number of Trotter steps the program is going to
  compute.

Optionally, you can specify the level of log messages for the program to report,
by setting `PHASE2_LOG` environment variable to be one of: `trace`,
`debug`, `info`, `warn`, `error`, `fatal`. E.g.,

```bash
mpirun -n 8 -x PHASE2_LOG=info ./ph2run-trott simul.h5 100
```

will compute 100 Trotter steps for a Hamiltonian specified in the file
`simul.h5` using 8 MPI processes, and write the result back to the same file.

See also the directory: [./simul](./simul) for an example of a simple
automated system.



Credits and License
-------------------

This software is distributed under the BSD-3Clause License. See [LICENSE](./LICENSE)
for more information.

Online repository available at: https://github.com/Quantum-for-Life/phase2


### Citation

Please cite this work as:

> Marek Miller, Jakob Gunther, Freek Witteveen, Matthew S. Teynor, Mihael Erakovic,
> Markus Reiher, Gemma C. Solomon, and Matthias Christandl,
> *phase2: Full-State Vector Simulation of Quantum Time Evolution at Scale*, arXiv preprint: TBA.


### Maintainers

* Marek Miller <mlm@math.ku.dk>

Report bugs or submit patches via [GitHub Issues].

[GitHub Issues]: https://github.com/Quantum-for-Life/phase2/issues
