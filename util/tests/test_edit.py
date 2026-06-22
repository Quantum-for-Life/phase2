"""ph2 strip, diff and attr (CLI specification)."""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")


@pytest.fixture
def solved(make_worksheet):
    """Worksheet with two sizeable results groups."""
    big = np.zeros((4096, 2))
    return make_worksheet(
        multidet=[(0, 1.0, 0.0)],
        results={"circ_trott": ({"delta": 0.1}, big),
                 "circ_qdrift": ({"step_size": 0.1, "depth": 4,
                                  "num_samples": 4096, "seed": 1},
                                 big)})


# -- strip ---------------------------------------------------------------

def test_strip_default(run_ph2, solved, tmp_path):
    size_before = solved.stat().st_size
    rc, _, err = run_ph2("strip", solved)
    assert rc == 0, err
    with h5py.File(solved) as f:
        assert "circ_trott" not in f and "circ_qdrift" not in f
        assert "pauli_hamil" in f and "state_prep" in f
    assert solved.stat().st_size < size_before
    assert not list(tmp_path.glob("*.tmp.*"))
    rc, out, _ = run_ph2("validate", solved)
    assert (rc, out) == (0, "")


def test_strip_preserves_payload_and_uuid(run_ph2, solved,
                                          tmp_path):
    with h5py.File(solved, "a") as f:
        f.attrs["uuid"] = "fixed-for-test"
    ref = tmp_path / "ref.h5"
    import shutil
    shutil.copy(solved, ref)
    run_ph2("strip", solved)
    with h5py.File(solved) as f:
        assert f.attrs["uuid"] == "fixed-for-test"
    rc, out, _ = run_ph2("diff", "--ignore-results", "--strict",
                         solved, ref)
    assert (rc, out) == (0, "")


def test_strip_noop_without_results(run_ph2, make_worksheet):
    ws = make_worksheet(multidet=[(0, 1.0, 0.0)])
    rc, _, _ = run_ph2("strip", ws)
    assert rc == 0


def test_strip_group_restricts(run_ph2, solved):
    rc, _, _ = run_ph2("strip", solved, "--group", "circ_trott")
    assert rc == 0
    with h5py.File(solved) as f:
        assert "circ_trott" not in f
        assert "circ_qdrift" in f


def test_strip_group_absent_exit1(run_ph2, solved):
    rc, _, err = run_ph2("strip", solved, "--group",
                         "circ_cmpsit")
    assert rc == 1
    assert "circ_cmpsit" in err


def test_strip_to_output(run_ph2, solved, tmp_path):
    out = tmp_path / "out.h5"
    rc, _, _ = run_ph2("strip", solved, "-o", out)
    assert rc == 0
    with h5py.File(solved) as f:
        assert "circ_trott" in f  # original untouched
    with h5py.File(out) as f:
        assert "circ_trott" not in f
    rc, _, err = run_ph2("strip", solved, "-o", out)
    assert rc == 2
    assert "exists" in err
    rc, _, _ = run_ph2("strip", solved, "-o", out, "--force")
    assert rc == 0


# -- diff ----------------------------------------------------------------

def test_diff_self(run_ph2, solved):
    rc, out, _ = run_ph2("diff", solved, solved)
    assert (rc, out) == (0, "")


def test_diff_tolerance(run_ph2, make_worksheet, tmp_path):
    import shutil
    a = make_worksheet(name="a.h5", multidet=[(0, 1.0, 0.0)])
    b = tmp_path / "b.h5"
    shutil.copy(a, b)
    with h5py.File(b, "a") as f:
        c = f["pauli_hamil/coeffs"][...]
        c[0] += 5e-13  # below the 1e-12 default tolerance
        del f["pauli_hamil/coeffs"]
        f["pauli_hamil"].create_dataset("coeffs", data=c)
    rc, out, _ = run_ph2("diff", a, b)
    assert (rc, out) == (0, "")
    rc, out, _ = run_ph2("diff", "--tol", "1e-15", a, b)
    assert rc == 1
    assert "dataset differs: /pauli_hamil/coeffs: max|d|=" in out


def test_diff_nan_padding_equal(run_ph2, solved, tmp_path):
    import shutil
    other = tmp_path / "o.h5"
    shutil.copy(solved, other)
    for p in (solved, other):
        with h5py.File(p, "a") as f:
            f["circ_trott/values"][-8:] = np.nan
    rc, out, _ = run_ph2("diff", solved, other)
    assert (rc, out) == (0, "")


def test_diff_attr_and_only_in(run_ph2, solved, tmp_path):
    import shutil
    other = tmp_path / "o.h5"
    shutil.copy(solved, other)
    with h5py.File(other, "a") as f:
        f["circ_trott"].attrs["delta"] = np.float64(0.2)
        del f["circ_qdrift"]
    rc, out, _ = run_ph2("diff", solved, other)
    assert rc == 1
    assert "attr differs: /circ_trott@delta: 0.1 != 0.2" in out
    assert "only in A: /circ_qdrift" in out
    rc, out, _ = run_ph2("diff", "--ignore-results", solved, other)
    assert (rc, out) == (0, "")


def test_diff_shape_and_dtype(run_ph2, make_worksheet, tmp_path):
    import shutil
    a = make_worksheet(name="a.h5")
    b = tmp_path / "b.h5"
    shutil.copy(a, b)
    with h5py.File(b, "a") as f:
        del f["pauli_hamil/coeffs"]
        f["pauli_hamil"].create_dataset(
            "coeffs", data=np.zeros(3, dtype=np.int64))
    rc, out, _ = run_ph2("diff", a, b)
    assert rc == 1
    assert "dtype differs: /pauli_hamil/coeffs" in out


def test_diff_uuid_ignored_unless_strict(run_ph2, tmp_path):
    terms = tmp_path / "T"
    terms.write_text("1.0 Z0\n")
    a, b = tmp_path / "a.h5", tmp_path / "b.h5"
    run_ph2("hamil", "paulis", terms, "-o", a)
    run_ph2("hamil", "paulis", terms, "-o", b)
    rc, out, _ = run_ph2("diff", a, b)
    assert (rc, out) == (0, "")
    rc, out, _ = run_ph2("diff", "--strict", a, b)
    assert rc == 1
    assert "/@uuid" in out


# -- attr ----------------------------------------------------------------

def test_attr_get_all_and_one(run_ph2, make_worksheet):
    ws = make_worksheet()
    rc, out, _ = run_ph2("attr", "get", ws, "/pauli_hamil")
    assert rc == 0
    assert "normalization=" in out and "offset=-1.0" in out
    rc, out, _ = run_ph2("attr", "get", ws, "/pauli_hamil",
                         "offset")
    assert (rc, out.strip()) == (0, "offset=-1.0")


def test_attr_get_missing_exit2(run_ph2, make_worksheet):
    ws = make_worksheet()
    rc, _, err = run_ph2("attr", "get", ws, "/pauli_hamil",
                         "nope")
    assert rc == 2
    assert "no attribute" in err
    rc, _, err = run_ph2("attr", "get", ws, "/absent")
    assert rc == 2


def test_attr_set_round_trip_and_dtype(run_ph2, make_worksheet):
    ws = make_worksheet()
    rc, _, _ = run_ph2("attr", "set", ws, "/pauli_hamil",
                       "offset", "-2.5")
    assert rc == 0
    rc, out, _ = run_ph2("attr", "get", ws, "/pauli_hamil",
                         "offset")
    assert out.strip() == "offset=-2.5"
    # dtype change f64 -> u64
    run_ph2("attr", "set", ws, "/pauli_hamil", "offset", "3",
            "--type", "u64")
    with h5py.File(ws) as f:
        assert f["pauli_hamil"].attrs["offset"].dtype == np.uint64


def test_attr_set_u8_keeps_validate_happy(run_ph2, tmp_path,
                                          testdata):
    import shutil
    ws = tmp_path / "w.h5"
    shutil.copy(testdata / "N4_closed.h5", ws)
    rc, _, _ = run_ph2("attr", "set", ws,
                       "/state_prep/coeff_matrix",
                       "closed_shell", "1", "--type", "u8")
    assert rc == 0
    rc, out, _ = run_ph2("validate", ws)
    assert (rc, out) == (0, "")


def test_attr_set_bad_value_exit2(run_ph2, make_worksheet):
    ws = make_worksheet()
    rc, _, err = run_ph2("attr", "set", ws, "/pauli_hamil",
                         "offset", "abc")
    assert rc == 2
    assert "cannot parse" in err
