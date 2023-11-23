# Specification: `simul.h5`

This document describes the format and the content of the _simulation
worksheet_ file. The file represents the state of the collective work of several
applications to obtain ground state energy estimation for a Hamiltonian,
modelling an electronic system in quantum chemistry by simulating
an algorithm related to quantum phase estimation (QPE).

# Format and file name

The format of the simulation file is HDF5.

The simulation file name is provided to the application as a command line
argument, as specified in the application help page or
documentation. If no such argument is given, the simulation file name
defaults to `simul.h5` in the working directory, i.e. the directory the
application was invoked from.

# Data types

For the description of data types, see the list
of [HDF5 Predefined Datatypes][hdf5-data-types].

# Structure

The simulation file contains named *groups* that consist of other
groups, *datasets*, i.e. homogenous multidimensional arrays of elements, and
*attributes*, i.e. metadata describing properties of dataset. The main group
of the simulation file is called the *root* group and has the name: `/`.

Names of other groups and attributes consist of only lowercase letters, numbers
and the underscore: `_`.

## Attribute: `/name`

- Type: string, `H5T_C_S1`
- Comment: Name given to the simulation run

## Attribute: `/circuit_name`

- Type: string, `H5T_C_S1`
- Comment: Name of a circuit (algorithm) simulated

## Group: `/fermion`

### Attribute: `/fermion/f2q_map`

- Type: string, `H5T_C_S1`
- Comment: Kind of fermion-to-qubit mapping used to obtain Pauli Hamiltonian.

### Dataset: `/fermion/fcidump`

[TBA]

## Group `/state_prep`

### Group: `/state_prep/multidet`

Contains a description of initial simulation state as a linear combination of
Slater determinants.

For given integers `NUM_TERMS, NUM_QUBITS >=1`:

- Dataset: `coeffs`
    - Type: double, `H5T_IEEE_F64LE`
    - Shape: `(NUN_TERMS,)`

- Dataset: `dets`
    - Type: unsigned char, `H5T_STD_U8LE`
    - Shape: `(NUM_TERMS, NUM_QUBITS)`
    - Comment: Rows denote computational basis states.

### Group: `/state_prep/mps`

Contains groups representing together a sequence of unitary
operations acting on specified parts of a quantum register. The sequence, when
applied in the correct order on the zero state of the quantum register, effects
the state preparation stage.

Each group must have a unique name starting with `unitary_`.

- Group: `/state_prep/mps/unitary_ID`
  (where `ID` is a unitary identifier, e.g. its index)

  For a given integer `N >= 2`:

    - Attribute: `index`:
        - Type: unsigned long, `H5T_STD_U64LE`
    - Dataset: `matrix`
        - Type: double, `H5T_IEEE_F64LE`
        - Shape: `(N,N)`
    - Dataset: `qubits`
        - Type: unsigned long, `H5T_STD_U64LE`
        - Shape: `(N,)`
        - Comment: Unique indices of qubit sites in the register to act on.

## Group: `/pauli_hamil`

For given integers `NUM_TERMS, NUM_QUBITS >=1`:

- Dataset: `coeffs`
    - Type: double, `H5T_IEEE_F64LE`
    - Shape: `(NUN_TERMS,)`


- Dataset: `paulis`
    - Type: unsigned char, `H5T_STD_U8LE`
    - Shape: `(NUM_TERMS, NUM_QUBITS)`
    - Comment: Elements in the dataset denote single-qubit Pauli operators
      according to the convention:
      ```text
      Id = 0, X = 1, Y = 2, Z = 3
      ```

## Group: `/time_series`

For a given integer `NUM_TERMS >= 1`:

- Dataset: `times`
    - Type: double, `H5T_IEEE_F64LE`
    - Shape: `(NUM_TERMS,)`


- Dataset: `values`
    - Type: double, `H5T_IEEE_F64LE`
    - Shape: `(NUM_TERMS, 2)`
    - Comment: Columns specify the real (column 1) and imaginary (column 2)
      part of a complex number. Uninitialized values (values to be computed)
      are designated by `NaN`.

## Group: `/log`

[TBA]

# Validation

[TBA]

[hdf5-data-types]: https://docs.hdfgroup.org/hdf5/v1_14/predefined_datatypes_tables.html