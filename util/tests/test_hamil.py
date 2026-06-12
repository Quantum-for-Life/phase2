"""ph2 hamil: TERMS grammar and FCIDUMP builder."""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

TERMS = """\
# H2-like toy
-0.4804   Z0
 0.3435   Z0 Z1   # trailing comment
 0.1809   X0 X1

 0.7137            # identity -> @offset
"""


@pytest.fixture
def terms(tmp_path):
    p = tmp_path / "TERMS"
    p.write_text(TERMS)
    return p


def _build(run_ph2, terms, out, *extra):
    return run_ph2("hamil", "paulis", terms, "-o", out, *extra)


def test_paulis_round_trip(run_ph2, terms, tmp_path):
    out = tmp_path / "h.h5"
    rc, _, _ = _build(run_ph2, terms, out)
    assert rc == 0
    with h5py.File(out) as f:
        g = f["pauli_hamil"]
        assert list(g["coeffs"][...]) == [-0.4804, 0.3435, 0.1809]
        assert g["paulis"][...].tolist() == [[3, 0], [3, 3], [1, 1]]
        s = abs(-0.4804) + 0.3435 + 0.1809
        assert g.attrs["normalization"] == pytest.approx(1 / s)
        assert g.attrs["offset"] == pytest.approx(0.7137)
        assert "uuid" in f.attrs
    rc, out_, _ = run_ph2("validate", "--hamil-only", out)
    assert (rc, out_) == (0, "")


def test_offset_flag_added(run_ph2, terms, tmp_path):
    out = tmp_path / "h.h5"
    _build(run_ph2, terms, out, "--offset", "-1.0")
    with h5py.File(out) as f:
        assert f["pauli_hamil"].attrs["offset"] == \
            pytest.approx(0.7137 - 1.0)


def test_identity_token_widens(run_ph2, tmp_path):
    p = tmp_path / "T"
    p.write_text("1.0 Z0 I3\n")
    out = tmp_path / "h.h5"
    rc, _, _ = _build(run_ph2, p, out)
    assert rc == 0
    with h5py.File(out) as f:
        assert f["pauli_hamil/paulis"].shape == (1, 4)


def test_duplicate_qubit_exit2(run_ph2, tmp_path):
    p = tmp_path / "T"
    p.write_text("1.0 Z0 X0\n")
    rc, _, err = _build(run_ph2, p, tmp_path / "h.h5")
    assert rc == 2
    assert "duplicate qubit" in err


def test_n_qubits_pad_and_floor(run_ph2, terms, tmp_path):
    out = tmp_path / "h.h5"
    _build(run_ph2, terms, out, "--n-qubits", "5")
    with h5py.File(out) as f:
        paulis = f["pauli_hamil/paulis"][...]
        assert paulis.shape == (3, 5)
        assert not paulis[:, 2:].any()
    rc, _, err = _build(run_ph2, terms, tmp_path / "i.h5",
                        "--n-qubits", "1")
    assert rc == 2
    assert "below the maximum" in err


def test_empty_hamiltonian_exit2(run_ph2, tmp_path):
    p = tmp_path / "T"
    p.write_text("# nothing\n0.5\n")
    rc, _, err = _build(run_ph2, p, tmp_path / "h.h5")
    assert rc == 2
    assert "empty Hamiltonian" in err


def test_bad_coefficient_exit2(run_ph2, tmp_path):
    p = tmp_path / "T"
    p.write_text("1+2j Z0\n")
    rc, _, err = _build(run_ph2, p, tmp_path / "h.h5")
    assert rc == 2
    assert "bad coefficient" in err


def test_sort_terms(run_ph2, tmp_path):
    p = tmp_path / "T"
    p.write_text("1.0 Z1\n2.0 X0\n")
    out = tmp_path / "h.h5"
    _build(run_ph2, p, out, "--sort-terms")
    with h5py.File(out) as f:
        # labels: Z1 -> "ZI", X0 -> "IX"; sorted: IX before ZI
        assert f["pauli_hamil/coeffs"][...].tolist() == [2.0, 1.0]


def test_overwrite_guard(run_ph2, terms, tmp_path):
    out = tmp_path / "h.h5"
    _build(run_ph2, terms, out)
    before = out.read_bytes()
    rc, _, err = _build(run_ph2, terms, out)
    assert rc == 2
    assert "exists" in err
    assert out.read_bytes() == before
    rc, _, _ = _build(run_ph2, terms, out, "--force")
    assert rc == 0
