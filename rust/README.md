# `phase2/rust`

To get Rust [nightly toolchain][rust-nightly]:

  ```bash
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  source ${HOME}/.cargo/env
  rustup toolchain install nightly
  rustup default nightly
  ```

[rust-nightly]: (https://rust-lang.github.io/rustup/concepts/channels.html)
