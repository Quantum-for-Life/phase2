# `simul` tree

This directory provides a basic automated system to perform the simulation 
with `ph2run`.

# Dependencies

To run the scripts, you'll need these python packages:

```bash
pip install h5py qiskit-nature matplotlib
```

# How to use it

1. Make sure `ph2run` compiled in `../build` 
(see [../README.md](../README.md) on how to compile the sources).
2. Put `FCIDUMP` file in this directory (with exactly this name).
3. Run:

```
make TIME_FACT=1.0 TROTT_STEPS=128 
```

or modify `Makefile` to your liking.