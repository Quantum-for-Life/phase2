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
documentation.

# Structure

The simulation file contains named *groups* that consist of other
groups, *datasets*, i.e. homogenous multidimensional arrays of elements, and
*attributes*, i.e. metadata describing properties of dataset. The main group
of the simulation file is called the *root* group and has the name: `/`.

Names of other groups and attributes consist of only lowercase letters, numbers
and the underscore: `_`.

## Group `/state_prep`

### Group: `/state_prep/multidet`

Contains a description of initial simulation state as a linear combination of
Slater determinants.

For given integers `NUM_TERMS, NUM_QUBITS >=1`:

- Dataset: `coeffs`
    - *Type*: `double`
    - *Shape*: `(NUN_TERMS,2)`
    - *Comment*: Columns specify the real (column 1) and imaginary (column 2)
      part of a complex number.


- Dataset: `dets`
    - *Type*: `unsigned char`
    - *Shape*: `(NUM_TERMS, NUM_QUBITS)`
    - *Comment*: Rows denote computational basis states.

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
        - *Type*: `unsigned long`
    - Dataset: `matrix`
        - *Type*: `double`
        - *Shape*: `(N,N)`
    - Dataset: `qubits`
        - *Type*: `unsigned long`
        - *Shape*: `(N,)`
        - *Comment*: Unique indices of qubit sites in the register to act on.

## Group: `/pauli_hamil`

For given integers `NUM_TERMS, NUM_QUBITS >=1`:

- Attribute: `offset`
    - *Type*: `double`
    - *Comment*: Coefficient multiplying the identity term in the
      Hamiltonian. Note that the term "`0 0 0 etc.`" must not appear in the
      dataset `paulis` below.

- Attribute: `normalization`
    - *Type*: `double`
    - *Comment*: Coefficients stored in dataset `coeffs` will be
      **multiplied** be this value and the simulation will run for the
      normalized Hamiltonian.

- Dataset: `coeffs`
    - *Type*: `double`
    - *Shape*: `(NUN_TERMS,)`

- Dataset: `paulis`
    - *Type*: `unsigned char`
    - *Shape*: `(NUM_TERMS, NUM_QUBITS)`
    - *Comment*: Elements in the dataset denote single-qubit Pauli operators
      according to the convention:

      ```text
      I = 0, X = 1, Y = 2, Z = 3
      ```

## Group: `/circ_trott`

- Attribute: `time_factor`
    - *Type*: `double`
    - *Comment*: Coefficient multiplying the time parameter for Hamiltonian
      simulation.

- Dataset: `values`
    - *Type*: `double`
    - *Shape*: `(NUM_STEPS, 2)`
    - *Comment*: Columns specify the real (column 1) and imaginary (column 2)
      part of a complex number. This dataset is created be the algorithm `silk`,
      and the value of `NUM_STEPS` is specified as a command-line argument.

## Group: `/circ_qdrift`

- Dataset: `values`
    - *Type*: `double`
    - *Shape*: `(NUM_STEPS, 2)`
    - *Comment*: Columns specify the real (column 1) and imaginary (column 2)
      part of a complex number. This dataset is created be the algorithm `silk`,
      and the value of `NUM_STEPS` is specified as a command-line argument.

[hdf5-data-types]: https://docs.hdfgroup.org/hdf5/v1_14/predefined_datatypes_tables.html

[uuid-rfc]: https://datatracker.ietf.org/doc/html/rfc4122
