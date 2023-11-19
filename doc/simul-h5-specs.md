# Specification: `simul.h5`

This document describes the format and the content of the _simulation
worksheet_ file. The file represents the state of the collective work of several
applications to obtain ground state energy estimation for a
Hamiltonian modelling a electronic system in quantum chemistry via
simulation of algorithms related to quantum phase estimation
(QPE).

# Definitions

The specification refers to the described files as the _simulation file_.

Any application working on the simulation file is in this document referred
to collectively as _the application_.

The simulation _data_ means the content of the simulation file. The data
can be:

- _well-defined_, i.e. structured correctly, according to this specification
- _unspecified_, i.e. having property that is beyond the scope of this
  specification. The application providing data that is unspecified
  here is required to document it explicitly.
- _undefined_, i.e. not conforming to this specification.

The application is free to assume that the data is well-defined at all times.
This allows optimisations and simplifies the program. Any undefined data must be
explicitly mentioned in by this specification.

*Example. It is an undefined situation when the length of the array
specifying the real coefficients of the Hamiltonian, the array of
Pauli strings to be multiplied by those coefficient do not match. Thus, the
application need not check if the arrays have the same
length and can always assume this property holds. It is free to
deduce the size of the Hamiltonian from either of the arrays without
performing bound checks. If the data was undefined this could lead to
incorrect result, a crash of the application or even corruption of the
data. It is therefore the responsibility of the application that
_writes_ into the simulation file to assure that the data is
well-defined.*

# Format and file name

The format of the simulation file is HDF5.

The simulation file name is provided to the application as a command line
argument, as specified in the application help page or
documentation. If no such argument is given, the simulation file name
defaults to `sumul.h5` in the working directory, i.e. the directory the
application was invoked from.

# Structure

The simulation file contains named *groups* that consist of other
groups, *datasets*, i.e. homogenous multidimensional arrays of elements, and
*attributes*, i.e. metadata describing properties of data. The main group
of the simulation file is called *root* group and has the name: `/`.

Names of other groups and attributes consist of only lowercase letters, numbers
and underscore: `_`.

## Attribute: `name`

## Attribute: `circuit`

## Group: `state_prep`

## Group: `pauli_hamil`

## Group: `time_series`

