# Specification: `simul.h5`

This document describes the format and the content of the _simulation
worksheet_ file. The file represents the state of the collective work of several
applications to obtain ground state energy estimation for a Hamiltonian,
modelling an electronic system in quantum chemistry by simulating
an algorithm related to quantum phase estimation (QPE).

The `/state_prep` group now supports two mutually exclusive
subtypes: the legacy flat multi-determinant list
(`/state_prep/multidet/`) and the dense Slater-Condon
coefficient matrix (`/state_prep/coeff_matrix/`).  Dispatch
between them happens at file open time via
`data_state_prep_kind`; see "Dispatch rules" below.

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

### Group: `/state_prep/coeff_matrix`

Carries a real coefficient matrix `C` (size `n_sites x n_occ`)
that the simulator expands into a dense superposition via the
Slater-Condon formula

    c(occ_alpha, occ_beta)
        = det(C_alpha[occ_alpha, :]) * det(C_beta[occ_beta, :])

at `circ_prepst()` time.

Attributes:

- `n_qubits` (`uint32`): total qubit count of the simulator
  register.  For an untapered state `n_qubits = 2 * n_sites`;
  for a tapered state `n_qubits = 2 * n_sites - 2`.
- `n_sites` (`uint32`): number of spatial orbitals `N`.
- `n_alpha` (`uint32`): occupation count of the alpha spin
  block; the column count of `C_alpha`.
- `n_beta`  (`uint32`): occupation count of the beta spin
  block; the column count of `C_beta`.
- `closed_shell` (`uint8`): `1` if `C_beta` is absent and
  `C_alpha` is reused for the beta block; `0` otherwise.
- `tapered` (`uint8`): `1` if bits 0 and `n_sites` must be
  dropped from each generated bitstring (Z₂ + S_z parity
  reduction); `0` otherwise.

Datasets:

- `C_alpha`: shape `(n_sites, n_alpha)`, dtype
  `H5T_IEEE_F64LE`.  Row `q` of the matrix gives the MO
  expansion coefficient for the alpha qubit at site `q`.
- `C_beta`: shape `(n_sites, n_beta)`, dtype
  `H5T_IEEE_F64LE`.  Present only if `closed_shell == 0`.

Optional subgroup `csf/` (CSF superposition):

- Attribute `n_components` (`uint32`): number of components.
  Must be `>= 1`; see the dispatch rules below.
- Subgroups `0/`, `1/`, ..., one per component.  Each
  carries:
    - Attribute `coefficient` (`double`): real weight of the
      component in the linear combination.
    - Datasets `C_alpha` and `C_beta` matching the top-level
      shape contract.

When `csf/` is present the simulator zeroes the register and
accumulates per-component contributions; the top-level
`C_alpha`/`C_beta` datasets are still required and are read
when the file is treated as a single-block state.

A `csf/` subgroup that exists but advertises `n_components
== 0` is rejected at load time (`data_coeff_matrix_load`
returns `-1` with a `log_error` line).  A single-block
state must omit the `csf/` subgroup outright rather than
declaring an empty superposition.

### Dispatch rules

The simulator probes both `/state_prep/multidet` and
`/state_prep/coeff_matrix` once at startup
(`data_state_prep_kind`).  Exactly one must be present:

| `/state_prep/multidet` | `/state_prep/coeff_matrix` | result |
|---|---|---|
| absent | absent | error: no state-prep subgroup found |
| present | absent | path: `STPREP_MULTIDET` |
| absent | present | path: `STPREP_COEFF_MATRIX` |
| present | present | error: ambiguous; rebuild `simul.h5` |

`STPREP_MULTIDET` and `STPREP_COEFF_MATRIX` are the enum
values from `enum stprep_kind` in `phase2/include/phase2/data.h`.
The both-present case is intentionally rejected at file-open
time so misbuilt paks never reach the time-evolution loop.

A `/state_prep/coeff_matrix/csf` subgroup that exists with
`n_components == 0` is also rejected (see the
`csf/` paragraph above): there is no silent fall-through to
the single-block path.  The three rejected configurations
for `coeff_matrix` are therefore: both subtypes present,
neither subtype present, and an explicit empty CSF
superposition.

### Worked example: N=4 closed-shell Huckel state

The same trial state can be expressed both as multidet and
coeff_matrix.  For `n_sites = 4`, `n_alpha = n_beta = 2`, the
matrix is `4 x 2` (8 real numbers); the multidet expansion
has up to `C(4,2) * C(4,2) = 36` non-zero amplitudes.
`phase2/test/data/N4_closed.h5` and
`phase2/test/data/N4_multidet.h5` carry the two encodings
of the same state; running one Trotter step on both
produces identical `/circ_trott/values` (see
`test/t-circ_trott_coeff.c`).

### Tapering

The tapered convention drops bits 0 and `n_sites` from every
generated bitstring before it is written into the register.
A tapered `multidet` file applies the same drop at write
time; the `coeff_matrix` path replays it at expand time per
generated bitstring, so the schema stays uniform with the
untapered case.  See `drop_two_bits` in
`phase2/state_prep_coeff.c` and the tapered fixture
`N8_tapered.h5`.

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

The simulator creates the group and pre-allocates the
`values` dataset with NaN-padded shape `(NUM_STEPS, 2)` at
init time; rank 0 then writes one row per Trotter step
plus `H5Fflush`.  A run that crashes or times out partway
leaves the file consistent: rows for completed steps carry
real numbers, the trailing rows remain NaN.  This applies
identically to `/circ_trott2`, `/circ_qdrift`, and
`/circ_cmpsit`.

- Attribute: `delta`
    - *Type*: `double`
  - *Comment*: Coefficient multiplying the time parameter in Hamiltonian
      simulation.

- Dataset: `values`
    - *Type*: `double`
    - *Shape*: `(NUM_STEPS, 2)`
    - *Comment*: Columns specify the real (column 1) and imaginary (column 2)
      part of a complex number. The value of `NUM_STEPS` is specified as a command-line argument.

## Group: `/circ_trott2`

Output of the symmetric (Strang) 2nd-order Trotter product
formula:

  S_2(δ) = ∏_{k=1..K} exp(-i h_k δ/2) · ∏_{k=K..1} exp(-i h_k δ/2)

Each step applies one forward sweep at `δ/2` followed by one
reverse sweep at `δ/2` over the same Hamiltonian terms.  See
`circ/trott2.c` for the implementation and `doc/phase2.md §5`
for the algorithm description.

- Attribute: `delta`
    - *Type*: `double`
    - *Comment*: Coefficient multiplying the time parameter
      in the symmetric Trotter step.

- Dataset: `values`
    - *Type*: `double`
    - *Shape*: `(NUM_STEPS, 2)`
    - *Comment*: Columns specify the real (column 1) and
      imaginary (column 2) part of a complex number.  The
      value of `NUM_STEPS` is specified as a command-line
      argument.

## Group: `/circ_qdrift`


- Attribute: `step_size`
    - *Type*: `double`
    - *Comment*: QDrift parameter: step size.

- Attribute: `num_samples`
    - *Type*: `unsigned long`
    - *Comment*: Number of independent samples in QDrift algorithm.

- Attribute: `depth`
    - *Type*: `unsigned long`
    - *Comment*: Circuit depth.

- Attribute: `seed`
    - *Type*: `unsigned long`
    - *Comment*: PRNG seed used for this run; non-zero.  Written by
      current `ph2run`; absent from files produced before it was
      recorded, so readers must treat it as optional.

- Dataset: `values`
    - *Type*: `double`
    - *Shape*: `(NUM_SAMPLES, 2)`
    - *Comment*: Columns specify the real (column 1) and imaginary (column 2)
      part of a complex number.

## Group: `/circ_cmpsit`

Output of the composite (partially randomised) 2nd-order
Trotter algorithm.  The Hamiltonian is split into a
deterministic top-`length` part applied lexicographically
and a randomised qDRIFT-style part of `depth` rotations per
step.  See `circ/cmpsit.c` and `doc/phase2.md §5` for the
algorithm description.

- Attribute: `length`
    - *Type*: `unsigned long`
    - *Comment*: Number of deterministic top-|c_k| terms.

- Attribute: `depth`
    - *Type*: `unsigned long`
    - *Comment*: Number of randomised terms drawn per step.

- Attribute: `angle_det`
    - *Type*: `double`
    - *Comment*: Step size for the deterministic part.

- Attribute: `angle_rand`
    - *Type*: `double`
    - *Comment*: Step size for the randomised part.

- Attribute: `steps`
    - *Type*: `unsigned long`
    - *Comment*: Number of composite Trotter steps per sample.

- Attribute: `seed`
    - *Type*: `unsigned long`
    - *Comment*: PRNG seed used for this run; non-zero.

- Dataset: `values`
    - *Type*: `double`
    - *Shape*: `(NUM_SAMPLES, 2)`
    - *Comment*: Columns specify the real (column 1) and
      imaginary (column 2) part of a complex number.

[hdf5-data-types]: https://docs.hdfgroup.org/hdf5/v1_14/predefined_datatypes_tables.html

[uuid-rfc]: https://datatracker.ietf.org/doc/html/rfc4122
