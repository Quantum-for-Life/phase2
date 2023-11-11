# qpe-simul-test

Test suite for [qpe-simul](https://github.com/Quantum-for-Life/qpe-simul). 🏭

# Installation

Tried under Linux only.

## Dependencies

- C/C++ toolchain and CMake:

  ```bash
  sudo apt install gcc g++ cmake
  ```

- OpenMPI:

  ```bash
  sudo apt install libopenmpi-dev openmpi-common
  ```

- Rust [nightly toolchain][rust-nightly]:

  ```bash
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  source ${HOME}/.cargo/env
  rustup toolchain install nightly
  rustup default nightly
  ```

[rust-nightly]: (https://rust-lang.github.io/rustup/concepts/channels.html)

## Getting the sources

Make sure you've got read access
to [github.com/Quantum-for-Life](https://github.com/Quantum-for-Life)
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

Then authenticate git:

```bash
gh auth login
```

Now you can clone the repo and its dependencies:

```bash
git clone --recurse-submodules https://github.com/Quantum-for-Life/qpe-simul-test.git
cd qpe-simul-test
```

Run the test suite by typing:

```bash
cargo test
```
