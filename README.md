# phase2

[![CI](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml)

Quantum phase estimation (QPE) with parallel computing.

# How to use it

Make sure you have `openmpi` installed. For example, on Ubuntu:

```bash
sudo apt install openmpi-common libopenmpi-dev
```

If you want to run the simulation on Euler cluster (ETHZ), load necessary
modules:

```bash
env2lmod
module load gcc
module load cmake
module load openmpi
```

Then

```bash
git clone --recurse-submodules git@github.com:Quantum-for-Life/qpe-simul.git
cd qpe-simul
mkdir build && cd build
cmake -DDISTRIBUTED=ON ..
make
```
