# ph2 — the worksheet toolkit

`util/ph2.py` is the command-line toolkit for phase2 worksheets
(`simul.h5`): inspect, validate, build, analyse and edit the HDF5
files specified in [simul-h5-specs.md](simul-h5-specs.md).  It runs
straight from the checkout — no installation step:

    util/ph2.py COMMAND [ARGS...]

## Dependencies

| Commands                       | Needs                                |
|--------------------------------|--------------------------------------|
| everything below               | `h5py`, `numpy` (`pip install -e ".[examples]"`) |
| `energy fft --peaks`           | + `scipy` (same extra)               |
| `hamil fcidump`                | + `qiskit-nature`, `pyscf` (`pip install -e ".[prep]"`) |

`ph2.py --help` works without any of them.

## Conventions

- Exit codes: **0** success, **1** semantic outcome (violations
  found, files differ, state-prep conflict), **2** usage / IO /
  missing dependency.
- Diagnostics go to stderr prefixed `ph2: `; machine-readable
  output (CSV, `name=value` lines) goes to stdout.
- Commands that write a new file (`-o`) refuse to overwrite an
  existing path unless `--force`.

---

## ph2 show — summarise a worksheet

    util/ph2.py show FILE...

One line per recognised top-level group; the `done/total` counter
of a results group counts the leading non-NaN rows of `values`
(completed steps or samples).

    $ util/ph2.py show examples/data/hamil.h5
    /pauli_hamil   251 terms, 10 qubits, norm=7.170948e-02, offset=-72.219274
    /state_prep    multidet: 1 det, 10 qubits

Unknown groups print `(unrecognised)`; legacy dangling soft links
print `(dangling link)`.  With several FILEs each block is headed
by `== FILE ==`.

## ph2 validate — check against the specification

    util/ph2.py validate [--hamil-only] FILE...

One violation per line, `FILE: H5PATH: message`; silent and exit 0
when clean.  Checks dataset dtypes and shapes, attribute presence,
Pauli bytes in {0,1,2,3} with no all-identity row, the exactly-one
state-prep rule, the CSF component contract,
`n_qubits == 2*n_sites - 2*tapered`, per-group result attributes
(`/circ_qdrift@seed` is optional), and that NaN rows form a suffix
of `values` (the crash-consistency guarantee).  `--hamil-only`
accepts intermediate paks without `/state_prep`.

    $ util/ph2.py validate test/data/N4_both.h5
    test/data/N4_both.h5: /state_prep: both multidet and coeff_matrix present

## ph2 hamil — build /pauli_hamil

    util/ph2.py hamil fcidump FCIDUMP -o OUT.h5 [--sort-terms] [--force]
    util/ph2.py hamil paulis  TERMS   -o OUT.h5 [--n-qubits N]
                              [--offset E0] [--sort-terms] [--force]

Both write a fresh file with a root `uuid` and `/pauli_hamil`:
`coeffs`, per-qubit `paulis` bytes (I=0, X=1, Y=2, Z=3),
`normalization` = 1/sum|c| (the simulator evolves a unit-norm
operator) and `offset` (the identity contribution, added back at
analysis time).  `--sort-terms` orders terms by Pauli label, which
maximises the simulator's batch-cache hit rate.

`fcidump` maps molecular integrals through qiskit-nature's
Jordan-Wigner transform; the identity term and the nuclear
repulsion fold into `offset`.

`paulis` reads a plain text file, one term per line; `#` starts a
comment, blank lines are ignored:

    # H2-like toy
    -0.4804   Z0
     0.3435   Z0 Z1
     0.1809   X0 X1
     0.7137            # bare coefficient = identity -> offset

Tokens are `[IXYZ]<qubit>`; a duplicate qubit within a term is an
error; `I<k>` only widens the implied register.  The qubit count
defaults to the highest index + 1 (`--n-qubits` may raise it).

    $ util/ph2.py hamil paulis TERMS -o toy.h5
    $ util/ph2.py validate --hamil-only toy.h5

## ph2 stprep — build /state_prep

    util/ph2.py stprep multidet INPUTST -f FILE [--no-reorder] [--force]
    util/ph2.py stprep coeff -f FILE --c-alpha CA [--c-beta CB]
                             [--tapered] [--csf W:CA[,CB]]... [--force]

Appends the trial state to an existing worksheet.  Both builders
refuse to create the invalid both-subtypes configuration (exit 1;
`--force` replaces the existing payload) and reject a payload whose
qubit count differs from `/pauli_hamil`.

`multidet` reads one determinant per line, `RE IM o0 o1 ...`, with
0/1 spin-orbital occupations interleaved alpha/beta (the
`simul.h5` convention stores alpha block then beta block; the
reorder happens here, `--no-reorder` for pre-blocked input):

    1.0 0.0  1 1 1 1 1 1 0 0 0 0

`coeff` writes a coefficient-matrix payload from CMAT files — one
row of `n_occ` floats per site, `#` comments and blank lines
ignored, all rows the same length:

    # (n_sites, n_occ) = (4, 2) alpha block
     0.6533  0.2706
     0.2706 -0.6533
    -0.2706 -0.6533
    -0.6533  0.2706

Dimensions are inferred from the shapes: `n_sites` = rows,
`n_alpha` = columns of CA, `n_beta` = columns of CB (closed shell —
no `--c-beta` — reuses CA).  `--tapered` marks Z2 + S_z tapering
(`n_qubits = 2*n_sites - 2`; the matrices are not transformed).
Each `--csf WEIGHT:CA[,CB]` adds a CSF component, numbered in flag
order; component shapes and shell must match the top level.

    $ util/ph2.py stprep coeff -f toy.h5 --c-alpha CA \
          --csf 0.6:CA --csf 0.8:CB

## ph2 energy — extract energies

    util/ph2.py energy fft  FILE [--group circ_trott|circ_trott2] [--peaks]
    util/ph2.py energy mc   FILE [--group circ_qdrift|circ_cmpsit]
    util/ph2.py energy ref  FILE
    util/ph2.py energy rpe  FILE [--group circ_trott|circ_trott2]
    util/ph2.py energy rpe-qdrift PREFIX --delta D --epsilon EPS

The single-file commands print one CSV line `E,E_ref,dE`: the
estimate, the trial-state energy `E_ref = <psi|H|psi>` computed
directly from the input (`nan` without `/state_prep`), and
`dE = E - E_ref`.  `ref` prints the bare `E_ref`; it covers both
state-prep subtypes (coeff_matrix via the Slater-Condon expansion,
CSF and tapering included — cost grows as `C(n_sites, n_alpha)^2`).

- **fft**: the Trotter series `values[s] = <psi|U^(s+1)|psi>` is a
  sum of tones; the dominant non-DC peak maps to
  `E = 2*pi*f / (delta * normalization) + offset`.  Resolution is
  one FFT bin = `2*pi/(steps*delta*norm)`; raise the step count to
  sharpen the line.  `--peaks` lists every detected `E,amplitude`.
- **mc**: the randomised drivers store independent samples at one
  effective time `T` (`depth*asin(step_size)` for qDRIFT,
  `steps*angle_det` for the composite);
  `E = arg(mean(values)) / (T * normalization) + offset`.  Trailing
  NaN rows are excluded; a decohered mean (|mean| < 0.1) warns.
- **rpe**: robust phase estimation on the dyadic subsequence
  `values[2^i - 1]`; `E0 = sqrt(1-x^2)/x * tan(x*theta) / norm /
  delta` with `x = 2^-J`.  Useful when the FFT resolution is not
  enough.  A NaN at a dyadic index aborts (run incomplete).
- **rpe-qdrift**: the same ladder over a depth sweep
  `PREFIX-0 .. PREFIX-J`, `J = ceil(log2(delta/epsilon))`, at fixed
  step size `x = 2^-J` (tan undoes the per-gate asin exactly).
  Prints `delta,epsilon,E0,E`.

```
$ util/ph2.py energy fft examples/simul/trott/simul.h5
-75.00018314380438,-74.96314677561794,-0.03703636818644
```

## ph2 strip — remove results groups

    util/ph2.py strip FILE [-o OUT.h5] [--group NAME]... [--force]

Shrinks a solved worksheet back to a runnable pak.  Default strips
every `circ_*` group present; `--group` restricts (a named group
that is absent is an error).  Removal is a copy-rewrite — h5py
deletion only unlinks and HDF5 does not return freed pages — either
to `-o OUT` or through a same-directory temp file that atomically
replaces FILE (this breaks hard links).  Root attributes survive;
dangling legacy links are dropped.

## ph2 diff — compare two worksheets

    util/ph2.py diff A B [--tol TOL] [--ignore-results] [--strict]

One difference per line; exit 0 equal, 1 differences, 2 error —
usable directly from Makefiles and CI:

    only in A: /circ_qdrift
    attr differs: /circ_trott@delta: 0.05 != 0.1
    shape differs: /pauli_hamil/coeffs: (251,) != (250,)
    dataset differs: /pauli_hamil/coeffs: max|d|=3.2e-09

Floats compare within `--tol` (absolute, default 1e-12) with
NaN == NaN, so identically NaN-padded `values` compare clean;
integers and strings compare exactly.  The root `uuid` is ignored
unless `--strict`; `--ignore-results` skips the `circ_*` groups.

## ph2 attr — read or write attributes

    util/ph2.py attr get FILE H5PATH [NAME]
    util/ph2.py attr set FILE H5PATH NAME VALUE
                         [--type f64|u64|u32|u8|str]

`get` prints `name=value` lines.  `set` writes with an explicit
dtype (default `f64`; the spec's unsigned-long attributes need
`u64`, the `closed_shell`/`tapered` flags `u8`); an existing
attribute is recreated so the dtype can change.  This is the sharp
tool — run `ph2 validate` afterwards:

    $ util/ph2.py attr set simul.h5 /circ_trott delta 0.05
    $ util/ph2.py validate simul.h5

---

## Worked example: FCIDUMP to energy

```sh
make                                          # build/ph2run/ph2run
util/ph2.py hamil fcidump FCIDUMP -o simul.h5 --sort-terms
util/ph2.py stprep multidet INPUTST -f simul.h5
util/ph2.py validate simul.h5
mpirun -n 2 build/ph2run/ph2run -S simul.h5 trott \
       --delta=0.1 --steps=4096
util/ph2.py energy fft simul.h5               # E,E_ref,dE
util/ph2.py strip simul.h5                    # back to a pak
```

The `examples/simul/*` Makefiles script exactly this flow over the
bundled water CAS(5,6) fixture; see
[../examples/TUTORIAL.md](../examples/TUTORIAL.md).
