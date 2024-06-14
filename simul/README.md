# `simul` tree

This directory provides a basic automated system to perform the simulation
with `ph2run`.

# Dependencies

To run the scripts, you'll need these python packages:

```bash
pip install h5py qiskit-nature pyscf
```

# How to use it

1. Make sure `ph2run-trott` is compiled in `../build`
   (see [../README.md](../README.md) on how to compile the sources).
2. Put `FCIDUMP` file in this directory (with exactly this name).
3. Put `INPUTST` file in the same directory.
4. Run:

```
make
```

or modify `Makefile` to your liking.

# `INPUTST` (input state) file format

*TODO: This section should be moved to `../doc` and merge with the rest of
the documentation.*

The format of the `INPUTST` file is pretty basic. Each line specifies one
state (Slater determinant) of a covex linear combination:

```text
F F I I I ... I
```

where `F` stands for a floating point number and `I` is either `0` or `1`.  
The two floating point numbers together denote real and imaginary part of a
complex coefficient in the linear combination. The integer numbers `I I ..
I` denote the occupation of electron orbitals. By convention, the *alpha*
and *gamma* type of orbitals interleave:

```latex
I_{1,\alpha} I_{1, \beta} I_{2,\alpha} I_{2, \beta} etc. 
```

The number of integer values `I I ... I` must be the same for all rows.

Whitespace and empty lines are ignored.
