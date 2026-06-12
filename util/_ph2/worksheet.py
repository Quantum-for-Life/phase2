"""Read-side worksheet summary (ph2 show).

One line per recognised top-level group; the "done" counter of a
results group is the number of leading non-NaN rows of its values
dataset (the spec guarantees NaN rows form a suffix, so this is the
number of completed steps/samples).
"""

import sys

import numpy as np

from . import RESULTS_GROUPS, open_ro

_W = 15  # name-column width


def _done(values):
    """Leading non-NaN row count of a (n, 2) values dataset."""
    nan = np.isnan(values[:, 0])
    idx = np.argmax(nan)
    if idx == 0 and not nan[0]:
        return len(values)
    return int(idx) if nan[idx] else len(values)


def _line_hamil(g):
    n, nqb = g["paulis"].shape
    return (f"{'/pauli_hamil':<{_W}}{n} terms, {nqb} qubits,"
            f" norm={g.attrs['normalization']:.6e},"
            f" offset={g.attrs['offset']:.6f}")


def _line_stprep(g):
    if "multidet" in g:
        m, nqb = g["multidet/dets"].shape
        det = "det" if m == 1 else "dets"
        body = f"multidet: {m} {det}, {nqb} qubits"
    elif "coeff_matrix" in g:
        cm = g["coeff_matrix"]
        shell = "closed" if cm.attrs["closed_shell"] else "open"
        body = (f"coeff_matrix: n_sites={cm.attrs['n_sites']}"
                f" n_alpha={cm.attrs['n_alpha']}"
                f" n_beta={cm.attrs['n_beta']} {shell}"
                f" tapered={cm.attrs['tapered']}")
        if "csf" in cm:
            ncomp = cm["csf"].attrs["n_components"]
            body += f", csf: {ncomp} components"
    else:
        body = "(no subtype)"
    return f"{'/state_prep':<{_W}}{body}"


def _line_results(name, g):
    done = _done(g["values"][...])
    total = g["values"].shape[0]
    a = g.attrs
    if name in ("circ_trott", "circ_trott2"):
        body = f"delta={a['delta']:g}, {done}/{total} steps"
    elif name == "circ_qdrift":
        body = f"step_size={a['step_size']:g} depth={a['depth']}"
        if "seed" in a:
            body += f" seed={a['seed']}"
        body += f", {done}/{total} samples"
    else:  # circ_cmpsit
        body = (f"length={a['length']} depth={a['depth']}"
                f" angle_det={a['angle_det']:g}"
                f" angle_rand={a['angle_rand']:g}"
                f" steps={a['steps']} seed={a['seed']},"
                f" {done}/{total} samples")
    return f"{'/' + name:<{_W}}{body}"


def summary(f):
    """Per-group summary lines for an open worksheet."""
    lines = []
    dangling = set()
    for name in f:
        try:
            f[name]
        except KeyError:  # legacy dangling soft link
            dangling.add(name)
    if "pauli_hamil" in f and "pauli_hamil" not in dangling:
        lines.append(_line_hamil(f["pauli_hamil"]))
    if "state_prep" in f and "state_prep" not in dangling:
        lines.append(_line_stprep(f["state_prep"]))
    for name in RESULTS_GROUPS:
        if name in f and name not in dangling:
            lines.append(_line_results(name, f[name]))
    known = {"pauli_hamil", "state_prep", *RESULTS_GROUPS}
    for name in f:
        if name in dangling:
            lines.append(f"{'/' + name:<{_W}}(dangling link)")
        elif name not in known:
            lines.append(f"{'/' + name:<{_W}}(unrecognised)")
    return lines


def cmd_show(args):
    rc = 0
    for i, path in enumerate(args.files):
        if len(args.files) > 1:
            if i:
                print()
            print(f"== {path} ==")
        try:
            with open_ro(path) as f:
                print("\n".join(summary(f)))
        except Exception as e:
            print(f"ph2: {e}", file=sys.stderr)
            rc = 2
    return rc
