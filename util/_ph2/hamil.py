"""/pauli_hamil builders (ph2 hamil).

Two sources, one output contract: a fresh HDF5 file carrying a root
`uuid` attribute and /pauli_hamil with per-term real `coeffs`,
per-qubit `paulis` bytes (I=0, X=1, Y=2, Z=3), and the scalar
attributes `normalization` (1/sum|c| over the non-identity terms, so
the simulator evolves a unit-norm operator) and `offset` (the
identity contribution, added back at analysis time).  The all-zero
Pauli row never enters the dataset (spec).

TERMS text grammar, one term per line ('#' starts a comment, blank
lines ignored):

    term  := COEFF token*
    COEFF := real float
    token := [IXYZ]<index>        e.g.  -0.4804  Z0 Z1

A bare coefficient is an identity term folded into `offset`; `I<k>`
tokens are accepted and only widen the implied qubit count;
duplicate qubit indices within a term are rejected.

The FCIDUMP path maps molecular integrals through qiskit-nature's
Jordan-Wigner transform (port of the retired
examples/scripts/parse_fcidump.py); qiskit-nature and pyscf are
imported lazily so the rest of ph2 works without them.
"""

import re
import uuid

import numpy as np

from . import Ph2Error, guard_out

_TOKEN = re.compile(r"^([IXYZ])([0-9]+)$")
_OPS = {"I": 0, "X": 1, "Y": 2, "Z": 3}


def parse_terms(lines, n_qubits=None, offset=0.0):
    """Parse TERMS lines -> (coeffs, ops, offset, n_qubits).

    ops is a list of {qubit: code} dicts, one per non-identity term.
    """
    coeffs, ops = [], []
    max_idx = -1
    for lineno, raw in enumerate(lines, 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        try:
            cf = float(fields[0])
        except ValueError:
            raise Ph2Error(
                f"TERMS line {lineno}: bad coefficient"
                f" {fields[0]!r}") from None
        term = {}
        for tok in fields[1:]:
            m = _TOKEN.match(tok)
            if not m:
                raise Ph2Error(f"TERMS line {lineno}:"
                               f" bad token {tok!r}")
            op, idx = m.group(1), int(m.group(2))
            if idx in term:
                raise Ph2Error(f"TERMS line {lineno}: duplicate"
                               f" qubit index {idx}")
            term[idx] = _OPS[op]
            max_idx = max(max_idx, idx)
        term = {q: c for q, c in term.items() if c}
        if term:
            coeffs.append(cf)
            ops.append(term)
        else:
            offset += cf
    if not ops:
        raise Ph2Error("TERMS: no non-identity terms"
                       " (empty Hamiltonian)")
    nqb = max_idx + 1
    if n_qubits is not None:
        if n_qubits < nqb:
            raise Ph2Error(f"--n-qubits {n_qubits} below the"
                           f" maximum used index ({max_idx})")
        nqb = n_qubits
    return coeffs, ops, offset, nqb


def _label(row):
    """Pauli label, highest qubit first (qiskit convention)."""
    return "".join("IXYZ"[b] for b in row[::-1])


def write_hamil(path, coeffs, paulis, offset, sort_terms, force):
    """Write a fresh worksheet with /pauli_hamil only."""
    import h5py
    guard_out(path, force)
    coeffs = np.asarray(coeffs, dtype=np.float64)
    paulis = np.asarray(paulis, dtype=np.uint8)
    if sort_terms:
        order = np.argsort([_label(r) for r in paulis])
        coeffs, paulis = coeffs[order], paulis[order]
    with h5py.File(path, "w") as f:
        f.attrs["uuid"] = str(uuid.uuid4())
        g = f.create_group("pauli_hamil")
        g.create_dataset("coeffs", data=coeffs)
        g.create_dataset("paulis", data=paulis)
        g.attrs["normalization"] = np.float64(
            1.0 / np.abs(coeffs).sum())
        g.attrs["offset"] = np.float64(offset)


def cmd_paulis(args):
    try:
        with open(args.terms) as fh:
            lines = fh.readlines()
    except OSError as e:
        raise Ph2Error(f"cannot read {args.terms}: {e}") from None
    coeffs, ops, offset, nqb = parse_terms(
        lines, n_qubits=args.n_qubits, offset=args.offset)
    paulis = np.zeros((len(ops), nqb), dtype=np.uint8)
    for r, term in enumerate(ops):
        for q, c in term.items():
            paulis[r, q] = c
    write_hamil(args.out, coeffs, paulis, offset,
                args.sort_terms, args.force)
    return 0


def cmd_fcidump(args):
    try:
        import qiskit_nature
        from qiskit_nature.second_q.formats.fcidump import FCIDump
        from qiskit_nature.second_q.formats.fcidump_translator \
            import fcidump_to_problem
        from qiskit_nature.second_q.mappers.jordan_wigner_mapper \
            import JordanWignerMapper
    except ImportError:
        raise Ph2Error(
            "hamil fcidump requires qiskit-nature and pyscf"
            " (pip install -e \".[prep]\")") from None

    # Suppress deprecation warnings in the FCIDump translation.
    qiskit_nature.settings.use_symmetry_reduced_integrals = True
    qiskit_nature.settings.use_pauli_sum_op = False

    problem = fcidump_to_problem(FCIDump.from_file(args.fcidump))
    op = JordanWignerMapper().map(problem.hamiltonian.second_q_op())

    # Qiskit labels are highest-qubit-first; reverse so column j is
    # qubit j.
    coeffs = np.array([c.real for c in op.coeffs], dtype=np.float64)
    paulis = np.array(
        [[_OPS[ch] for ch in p.to_label()[::-1]] for p in op.paulis],
        dtype=np.uint8)

    # Fold the identity term and the nuclear repulsion into offset;
    # the all-zero row must not enter the dataset.
    nuc = problem.nuclear_repulsion_energy or 0.0
    ident = ~paulis.any(axis=1)
    offset = float(coeffs[ident].sum()) + nuc
    write_hamil(args.out, coeffs[~ident], paulis[~ident], offset,
                args.sort_terms, args.force)
    return 0
