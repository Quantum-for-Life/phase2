# phase2

[![CI](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Quantum-for-Life/phase2/actions/workflows/CI.yml)

Quantum phase estimation (QPE) with parallel computing. 🏭

# Dependencies

To run the simulations, you will need:

- Linux platform, with a C/C++ compiler toolchain
- CMake>=3.9
- [HDF5][hdf5-website]
- optional: [OpenMPI][openmpi-website]

To prepare input for the simulation and to run the test suite, we require:

- Python >= 3.10 and a virtual environment with several packages installed
  locally (see [README.md](./python/README.md) in the `python/` subsystem for
  more details.
- Rust [nightly toolchain][rust-nightly] (see [README.md](./rust/README.md) in
  `rust/` for more info).

Assuming we're on Ubuntu, to install the dependencies, run:

```bash
sudo apt install gcc g++ cmake 
sudo apt install libhdf5-dev hdf5-tools

sudo apt install python3 python3-dev python3-venv

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ${HOME}/.cargo/env
rustup toolchain install nightly
rustup default nightly
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
module load openmpi
```

[hdf5-website]: https://www.hdfgroup.org/solutions/hdf5/

[openmpi-website]: https://www.open-mpi.org/

[rust-nightly]: (https://rust-lang.github.io/rustup/concepts/channels.html)

## Getting the sources

Make sure you've got read access
to [https://github.com/Quantum-for-Life](https://github.com/Quantum-for-Life)
and you can clone repositories via HTTPS. The easiest way to do it is to
first install [GitHub CLI](https://cli.github.com/):

```bash
type -p curl >/dev/null || (sudo apt update && sudo apt install curl -y)
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | sudo dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg \
&& sudo chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
&& echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | sudo tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
&& sudo apt update \
&& sudo apt install gh -y
```

Then authenticate Git:

```bash
gh auth login
```

Now you can clone the repo and its dependencies:

```bash
git clone --recurse-submodules https://github.com/Quantum-for-Life/phase2.git
```

If you prefer to access GitHub via SSH instead, be aware the repositories
this package depends on as submodules have their URLs hard-coded as HTTPS.
And some of them are still private to _Quantum-for-Life_. You will need to
change the remote repo URL for each submodule to its SSH equivalent.

# Compiling the source code

```bash
cd phase2
mkdir build && cd build
cmake ..
make
```

If you want to enable MPI support, run CMake with `DISTRIBUTED` flag set:

```bash
cmake -DDISTRIBUTED=ON ..
make
```

# Testing

Run the test suite by typing:

```bash
cargo test
```

# How to use it

[TBA]
