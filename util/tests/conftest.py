"""Shared fixtures for the ph2 test suite.

Makes the private _ph2 package importable, exposes the repository
fixture paths, provides builders for synthetic worksheets in
tmp_path, and a subprocess runner for CLI-level tests.
"""

import subprocess
import sys
from pathlib import Path

import pytest

UTIL = Path(__file__).resolve().parents[1]
REPO = UTIL.parent
sys.path.insert(0, str(UTIL))


@pytest.fixture
def repo():
    return REPO


@pytest.fixture
def testdata():
    return REPO / "test" / "data"


@pytest.fixture
def examples_data():
    return REPO / "examples" / "data"


@pytest.fixture
def run_ph2():
    """Run util/ph2.py as a subprocess: (rc, stdout, stderr)."""
    def run(*args):
        r = subprocess.run(
            [sys.executable, str(UTIL / "ph2.py"), *map(str, args)],
            capture_output=True, text=True)
        return r.returncode, r.stdout, r.stderr
    return run


@pytest.fixture
def make_worksheet(tmp_path):
    """Build a synthetic worksheet; returns its path.

    Keyword groups: hamil (default on), multidet (list of
    (index, re, im)), results (dict group -> (attrs, values)).
    """
    h5py = pytest.importorskip("h5py")
    import numpy as np

    def make(name="ws.h5", hamil=True, nqb=2, multidet=None,
             results=None):
        path = tmp_path / name
        with h5py.File(path, "w") as f:
            if hamil:
                g = f.create_group("pauli_hamil")
                paulis = np.zeros((2, nqb), dtype=np.uint8)
                paulis[0, 0] = 3          # Z0
                paulis[1, :2] = (1, 1)    # X0 X1
                g.create_dataset("coeffs",
                                 data=np.array([0.5, 0.25]))
                g.create_dataset("paulis", data=paulis)
                g.attrs["normalization"] = np.float64(1.0 / 0.75)
                g.attrs["offset"] = np.float64(-1.0)
            if multidet is not None:
                sp = f.require_group("state_prep")
                md = sp.create_group("multidet")
                m = len(multidet)
                coeffs = np.zeros((m, 2))
                dets = np.zeros((m, nqb), dtype=np.uint8)
                for r, (idx, re, im) in enumerate(multidet):
                    coeffs[r] = (re, im)
                    for j in range(nqb):
                        dets[r, j] = (idx >> j) & 1
                md.create_dataset("coeffs", data=coeffs)
                md.create_dataset("dets", data=dets)
            for gname, (attrs, values) in (results or {}).items():
                g = f.create_group(gname)
                for k, v in attrs.items():
                    g.attrs[k] = v
                g.create_dataset("values", data=np.asarray(values))
        return path

    return make
