#!/usr/bin/env python3
# t-ref-coeff_matrix.py
#
# Strict cross-validation harness.  For every shipped
# coeff_matrix fixture:
#
#   1. Invoke the C-side Slater-Condon expansion via
#      `t-state_prep_coeff_expand --dump <out> <fixture>`
#      and parse the resulting (idx, re, im) lines.
#   2. Run the independent Python reference at
#      `ref/coeff_matrix_reference.py` against the same C
#      matrices read straight from the fixture H5.
#   3. Compare the two (idx, amplitude) dicts term by
#      term: identical index sets and per-index abs-diff
#      <= TOL.
#
# Exits 0 iff every fixture matches.  On mismatch, prints
# the fixture name, the worst-case index, both
# amplitudes, and the diff, then exits 1.

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "ref"))
from coeff_matrix_reference import expand_coeff_matrix  # noqa: E402

TOL = 1e-12

FIXTURES = [
    "N4_closed",
    "N4_open",
    "N4_csf",
    "N8_untapered",
    "N8_tapered",
]


def load_C(path: Path):
    """Read the full coeff_matrix description from an H5 fixture.

    Returns ``(attrs, top_Ca, top_Cb, csf_blocks)`` where
    ``csf_blocks`` is ``None`` for single-block fixtures and a
    list ``[(coefficient, Ca_k, Cb_k), ...]`` for CSF fixtures.
    """
    with h5py.File(path, "r") as f:
        g = f["state_prep/coeff_matrix"]
        attrs = dict(g.attrs)
        Ca = np.asarray(g["C_alpha"][...], dtype=np.float64)
        Cb = (
            np.asarray(g["C_beta"][...], dtype=np.float64)
            if "C_beta" in g
            else None
        )
        csf_blocks = None
        if "csf" in g:
            csf = g["csf"]
            n_components = int(csf.attrs["n_components"])
            csf_blocks = []
            for k in range(n_components):
                sub = csf[str(k)]
                cf = float(sub.attrs["coefficient"])
                Ca_k = np.asarray(sub["C_alpha"][...],
                                  dtype=np.float64)
                Cb_k = (
                    np.asarray(sub["C_beta"][...],
                               dtype=np.float64)
                    if "C_beta" in sub else None
                )
                csf_blocks.append((cf, Ca_k, Cb_k))
        return attrs, Ca, Cb, csf_blocks


def c_dump(binary: Path, fixture: Path) -> dict[int, complex]:
    """Invoke the C dump binary and parse its output."""
    with tempfile.NamedTemporaryFile(
        mode="w+", suffix=".amps", delete=False
    ) as tf:
        tmp_path = Path(tf.name)
    try:
        rc = subprocess.run(
            [str(binary), "--dump", str(tmp_path), str(fixture)],
            check=False,
        )
        if rc.returncode != 0:
            raise RuntimeError(
                f"{binary} --dump {fixture} returned {rc.returncode}"
            )
        out: dict[int, complex] = {}
        with tmp_path.open() as fh:
            for line in fh:
                parts = line.split()
                if len(parts) != 3:
                    raise RuntimeError(f"malformed dump line: {line!r}")
                idx = int(parts[0])
                re = float(parts[1])
                im = float(parts[2])
                out[idx] = complex(re, im)
        return out
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def py_expand(fixture: Path) -> dict[int, complex]:
    """Run the vendored reference against the fixture C matrices.

    Mirrors phase2's ``state_prep_coeff_expand_all`` dispatch: when
    a ``csf/`` subgroup is present, the result is the coefficient-
    weighted sum over per-component expansions; otherwise it is the
    single-block expansion of the top-level matrices.
    """
    attrs, Ca, Cb, csf_blocks = load_C(fixture)
    n_sites = int(attrs["n_sites"])
    n_alpha = int(attrs["n_alpha"])
    n_beta  = int(attrs["n_beta"])
    tapered = bool(int(attrs["tapered"]))
    if csf_blocks is None:
        real = expand_coeff_matrix(
            n_sites, n_alpha, n_beta, Ca, Cb, tapered,
        )
    else:
        real = {}
        for cf, Ca_k, Cb_k in csf_blocks:
            block = expand_coeff_matrix(
                n_sites, n_alpha, n_beta, Ca_k, Cb_k, tapered,
                weight=cf,
            )
            for idx, amp in block.items():
                real[idx] = real.get(idx, 0.0) + amp
        real = {i: a for i, a in real.items() if abs(a) >= 1e-12}
    return {idx: complex(amp, 0.0) for idx, amp in real.items()}


def compare(name: str, c_amps, py_amps) -> bool:
    """Return True if (name) passed.  Prints diagnostics on failure."""
    c_keys = set(c_amps)
    py_keys = set(py_amps)
    only_c = c_keys - py_keys
    only_py = py_keys - c_keys
    if only_c or only_py:
        print(
            f"FAIL[{name}]: index set differs;"
            f" C-only={sorted(only_c)[:5]}"
            f" Py-only={sorted(only_py)[:5]}"
            f" (sizes: C={len(c_keys)}, Py={len(py_keys)})",
            file=sys.stderr,
        )
        return False
    worst_idx = -1
    worst_diff = 0.0
    for idx in c_keys:
        d = abs(c_amps[idx] - py_amps[idx])
        if d > worst_diff:
            worst_diff = d
            worst_idx = idx
    if worst_diff > TOL:
        print(
            f"FAIL[{name}]: worst-case idx={worst_idx}"
            f" C={c_amps[worst_idx]} Py={py_amps[worst_idx]}"
            f" diff={worst_diff:.3e} > {TOL:.0e}",
            file=sys.stderr,
        )
        return False
    print(
        f"{name}: {len(c_keys)} amps match"
        f" (worst diff={worst_diff:.2e})"
    )
    return True


def main():
    # The build places compiled artefacts under build/<srcdir>/.
    # HERE is test/, so the helper binary lives at
    # <repo>/build/test/t-state_prep_coeff_expand.
    binary = HERE.parent / "build" / "test" / "t-state_prep_coeff_expand"
    if not binary.is_file() or not os.access(binary, os.X_OK):
        print(
            f"missing or non-executable {binary};"
            " run `make build-test` first",
            file=sys.stderr,
        )
        return 2

    base = HERE / "data"
    failed = []
    for name in FIXTURES:
        fixture = base / f"{name}.h5"
        if not fixture.exists():
            print(f"missing fixture {fixture}", file=sys.stderr)
            failed.append(name)
            continue
        c_amps = c_dump(binary, fixture)
        py_amps = py_expand(fixture)
        if not compare(name, c_amps, py_amps):
            failed.append(name)

    if failed:
        print(
            f"FAIL: {len(failed)} fixture(s) mismatched:"
            f" {failed}",
            file=sys.stderr,
        )
        return 1
    print("all fixtures: strict cross-validation OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
