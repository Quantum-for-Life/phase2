"""/state_prep builders and the expansion oracle (ph2 stprep).

Two builders append the state-prep payload to an existing worksheet,
with the exact dtypes of the committed fixtures (uint32 counts,
uint8 flags, float64 data):

  multidet      from an INPUTST text file: one determinant per line,
                `RE IM o0 o1 ...` with 0/1 occupations interleaved
                alpha/beta, reordered to alpha-block | beta-block.
  coeff_matrix  from CMAT text files (one row of n_occ floats per
                site, '#' comments); dimensions are inferred from
                the matrix shapes, closed shell from the absence of
                a beta block, and CSF components from repeated
                --csf WEIGHT:CA[,CB] flags.

Both refuse to create the spec's both-present configuration: a
worksheet whose /state_prep already carries a subtype exits 1
unless --force replaces it (copy-rewrite).  If /pauli_hamil exists,
its qubit width must match the payload.

expand_coeff_matrix is the Slater-Condon expansion oracle, ported
verbatim from test/ref/coeff_matrix_reference.py (this package must
not import from test/); a parity test pins the port to the
original.  Bitstring convention: alpha qubits at bits [0, n_sites),
beta at [n_sites, 2*n_sites); tapering drops bit 0 and bit n_sites.
"""

from itertools import combinations

import numpy as np

from . import Ph2Error, Ph2Failure, copy_without

PRUNE = 1e-12


def expand_coeff_matrix(n_sites, n_alpha, n_beta, C_alpha, C_beta,
                        tapered, weight=1.0):
    """Slater-Condon expansion -> {basis index: amplitude}."""
    Cb = C_beta if C_beta is not None else C_alpha
    out = {}

    alpha = []
    for occ in combinations(range(n_sites), n_alpha):
        d = 1.0 if n_alpha == 0 else \
            float(np.linalg.det(C_alpha[list(occ), :]))
        alpha.append((occ, d))

    beta = []
    for occ in combinations(range(n_sites), n_beta):
        d = 1.0 if n_beta == 0 else \
            float(np.linalg.det(Cb[list(occ), :]))
        beta.append((occ, d))

    for occ_a, da in alpha:
        for occ_b, db in beta:
            cf = weight * da * db
            if abs(cf) < PRUNE:
                continue
            idx = 0
            for q in occ_a:
                idx |= 1 << q
            for q in occ_b:
                idx |= 1 << (n_sites + q)
            if tapered:
                low = (idx >> 1) & ((1 << (n_sites - 1)) - 1)
                high = idx >> (n_sites + 1)
                idx = low | (high << (n_sites - 1))
            out[idx] = out.get(idx, 0.0) + cf
    return out


def load_cmat(path):
    """Parse a CMAT file -> float64 (n_sites, n_occ) array."""
    rows = []
    try:
        with open(path) as fh:
            for lineno, raw in enumerate(fh, 1):
                line = raw.split("#", 1)[0].strip()
                if not line:
                    continue
                try:
                    rows.append([float(x) for x in line.split()])
                except ValueError:
                    raise Ph2Error(f"{path} line {lineno}:"
                                   " bad matrix entry") from None
    except OSError as e:
        raise Ph2Error(f"cannot read {path}: {e}") from None
    if not rows:
        raise Ph2Error(f"{path}: empty matrix")
    if len({len(r) for r in rows}) != 1:
        raise Ph2Error(f"{path}: ragged rows")
    return np.array(rows, dtype=np.float64)


def parse_inputst(path):
    """Parse INPUTST -> (coeffs (M,2), occupations (M, nqb))."""
    coeffs, dets = [], []
    try:
        with open(path) as fh:
            for lineno, raw in enumerate(fh, 1):
                fields = raw.split()
                if not fields:
                    continue
                if len(fields) < 3:
                    raise Ph2Error(f"{path} line {lineno}:"
                                   " expected RE IM o0 o1 ...")
                try:
                    re_, im = float(fields[0]), float(fields[1])
                    occ = [int(x) for x in fields[2:]]
                except ValueError:
                    raise Ph2Error(f"{path} line {lineno}:"
                                   " bad field") from None
                if any(o not in (0, 1) for o in occ):
                    raise Ph2Error(f"{path} line {lineno}:"
                                   " occupation not 0/1")
                coeffs.append((re_, im))
                dets.append(occ)
    except OSError as e:
        raise Ph2Error(f"cannot read {path}: {e}") from None
    if not dets:
        raise Ph2Error(f"{path}: no determinants")
    if len({len(d) for d in dets}) != 1:
        raise Ph2Error(f"{path}: occupation counts differ"
                       " across rows")
    return (np.array(coeffs, dtype=np.float64),
            np.array(dets, dtype=np.uint8))


def _reorder(dets):
    """Interleaved alpha/beta columns -> alpha block | beta block."""
    nqb = dets.shape[1]
    if nqb % 2:
        raise Ph2Error("interleaved INPUTST needs an even"
                       " occupation count (use --no-reorder?)")
    order = list(range(0, nqb, 2)) + list(range(1, nqb, 2))
    return dets[:, order]


def _open_for_stprep(path, nqb, force):
    """Open FILE for appending /state_prep; enforce the guards."""
    import h5py
    try:
        f = h5py.File(path, "a")
    except OSError as e:
        raise Ph2Error(f"cannot open {path}: {e}") from None
    if "state_prep" in f:
        sp = f["state_prep"]
        if "multidet" in sp or "coeff_matrix" in sp:
            if not force:
                f.close()
                raise Ph2Failure(
                    f"{path}: /state_prep already carries a"
                    " subtype (use --force to replace)")
            f.close()
            copy_without(path, {"state_prep"})
            f = h5py.File(path, "a")
    if "pauli_hamil" in f:
        width = f["pauli_hamil/paulis"].shape[1]
        if width != nqb:
            f.close()
            raise Ph2Failure(
                f"{path}: /pauli_hamil has {width} qubits,"
                f" state prep has {nqb}")
    return f


def cmd_multidet(args):
    coeffs, dets = parse_inputst(args.inputst)
    if not args.no_reorder:
        dets = _reorder(dets)
    nqb = dets.shape[1]
    with _open_for_stprep(args.file, nqb, args.force) as f:
        md = f.require_group("state_prep").create_group("multidet")
        md.create_dataset("coeffs", data=coeffs)
        md.create_dataset("dets", data=dets)
    return 0


def _parse_csf(spec, closed):
    """--csf WEIGHT:CA[,CB] -> (weight, C_alpha, C_beta_or_None)."""
    try:
        weight, paths = spec.split(":", 1)
        weight = float(weight)
    except ValueError:
        raise Ph2Error(f"--csf {spec!r}: expected"
                       " WEIGHT:CA[,CB]") from None
    parts = paths.split(",")
    if closed and len(parts) != 1 or not closed and len(parts) != 2:
        raise Ph2Error(f"--csf {spec!r}: component shell differs"
                       " from the top level")
    ca = load_cmat(parts[0])
    cb = load_cmat(parts[1]) if len(parts) == 2 else None
    return weight, ca, cb


def cmd_coeff(args):
    Ca = load_cmat(args.c_alpha)
    n_sites, n_alpha = Ca.shape
    closed = args.c_beta is None
    if closed:
        Cb, n_beta = None, n_alpha
    else:
        Cb = load_cmat(args.c_beta)
        if Cb.shape[0] != n_sites:
            raise Ph2Error(f"{args.c_beta}: row count differs"
                           f" from {args.c_alpha}")
        n_beta = Cb.shape[1]
    csf = []
    for spec in args.csf or ():
        weight, ca_k, cb_k = _parse_csf(spec, closed)
        if ca_k.shape != Ca.shape or \
                (cb_k is not None and cb_k.shape != Cb.shape):
            raise Ph2Error(f"--csf {spec!r}: component shape"
                           " differs from the top level")
        csf.append((weight, ca_k, cb_k))
    tapered = 1 if args.tapered else 0
    nqb = 2 * n_sites - 2 * tapered

    with _open_for_stprep(args.file, nqb, args.force) as f:
        sp = f.require_group("state_prep")
        cm = sp.create_group("coeff_matrix")
        cm.attrs["n_qubits"] = np.uint32(nqb)
        cm.attrs["n_sites"] = np.uint32(n_sites)
        cm.attrs["n_alpha"] = np.uint32(n_alpha)
        cm.attrs["n_beta"] = np.uint32(n_beta)
        cm.attrs["closed_shell"] = np.uint8(1 if closed else 0)
        cm.attrs["tapered"] = np.uint8(tapered)
        cm.create_dataset("C_alpha", data=Ca)
        if not closed:
            cm.create_dataset("C_beta", data=Cb)
        if csf:
            g = cm.create_group("csf")
            g.attrs["n_components"] = np.uint32(len(csf))
            for k, (weight, ca_k, cb_k) in enumerate(csf):
                blk = g.create_group(str(k))
                blk.attrs["coefficient"] = np.float64(weight)
                blk.create_dataset("C_alpha", data=ca_k)
                if not closed:
                    blk.create_dataset("C_beta", data=cb_k)
    return 0
