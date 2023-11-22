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
defaults to `sumul.h5` in the working directory, i.e. the directory the
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

## Group: `/state_prep`

This group contains datasets representing a sequence of unitary operations.
The sequence, when applied in the correct order on
the zero state of a quantum register, effects the state preparation stage.

Each data set in the group has an attribute (`index = 0,1,2,..`), which
specifies the order of applying the operator on the input state.

- Dataset
    - Type: double, `H5T_IEEE_F64LE`
    - Shape: `(N,N)`, where `N>= 2`
    - Attribute: `index`:
        - Type: unsigned long, `H5T_STD_U64LE`

## Group: `/pauli_hamil`

For given integers `NUM_TERMS, NUM_QUBITS >=1`:

- Dataset: `coeffs`
    - Type: double, `H5T_IEEE_F64LE`
    - Shape: `(NUN_TERMS, 1)`


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
    - Shape: `(NUM_TERMS, 1)`


- Dataset: `values`
    - Type: double, `H5T_IEEE_F64LE`
    - Shape: `(NUM_TERMS, 2)`
    - Comment: Columns specify the real (column 1) and imaginary (column 2)
      part of a complex number. Uninitialized values (values to be computed)
      are designated by `NaN`.

[hdf5-data-types]: https://docs.hdfgroup.org/hdf5/v1_14/predefined_datatypes_tables.html