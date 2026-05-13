# `simul` tree

This directory provides a basic automated system to perform the simulation
with `ph2run`.

# Dependencies

To run the scripts, you'll need these python packages:

```bash
pip install h5py qiskit-nature pyscf
```

# How to use it

1. Make sure `ph2run` is compiled at `../ph2run/ph2run`
   (see [../README.md](../README.md) on how to compile the
   sources).
2. Put an `FCIDUMP` file in this directory (with exactly
   this name).
3. Put an `INPUTST` file in the same directory.
4. Run:

```
make
```

or modify `Makefile` to your liking.

# `INPUTST` (input state) file format

The format of the `INPUTST` file is plain text.  Each line
specifies one Slater determinant of a convex linear
combination:

```text
F F I I I ... I
```

where `F` is a floating-point number and `I` is either `0`
or `1`.  The two floats give the real and imaginary part of
the complex coefficient in the linear combination.  The
integers `I I ... I` denote occupation of spin orbitals.
By convention, the *alpha* and *beta* type orbitals
interleave:

```latex
I_{1,\alpha} I_{1, \beta} I_{2,\alpha} I_{2, \beta} etc.
```

The number of integer values `I I ... I` must be the same
for all rows.  Whitespace and empty lines are ignored.
