"""ph2 validate: one test per specification check."""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

GOOD = ["N4_closed.h5", "N4_multidet.h5", "N8_tapered.h5",
        "N8_untapered.h5"]


@pytest.mark.parametrize("name", GOOD)
def test_clean_fixture(run_ph2, testdata, name):
    rc, out, _ = run_ph2("validate", testdata / name)
    assert (rc, out) == (0, "")


def test_clean_examples_hamil(run_ph2, examples_data):
    rc, out, _ = run_ph2("validate", examples_data / "hamil.h5")
    assert (rc, out) == (0, "")


@pytest.mark.parametrize("name", ["N4_open.h5", "N4_csf.h5"])
def test_stprep_only_fixture(run_ph2, testdata, name):
    # These fixtures intentionally ship without a Hamiltonian; the
    # state-prep payload itself must be the only clean part.
    rc, out, _ = run_ph2("validate", testdata / name)
    assert rc == 1
    lines = out.strip().splitlines()
    assert len(lines) == 1
    assert "/pauli_hamil: missing" in lines[0]


def test_legacy_dangling_links(run_ph2, testdata):
    # H2O_CAS56.h5 carries legacy soft links circ_trott{,er} ->
    # "trotter_steps" that resolve to nothing; ph2run unlinks such
    # junk on group create, validate reports it.
    rc, out, _ = run_ph2("validate", testdata / "H2O_CAS56.h5")
    assert rc == 1
    assert "/circ_trott: dangling link" in out
    assert "/circ_trotter: dangling link" in out


def test_both_subtypes(run_ph2, testdata):
    rc, out, _ = run_ph2("validate", testdata / "N4_both.h5")
    assert rc == 1
    assert "both multidet and coeff_matrix" in out


def test_csf_empty(run_ph2, testdata):
    rc, out, _ = run_ph2("validate", testdata / "N4_csf_empty.h5")
    assert rc == 1
    assert "n_components" in out


# -- synthetic negatives ------------------------------------------------

def _edit(path):
    return h5py.File(path, "a")


@pytest.fixture
def ws(make_worksheet):
    """Valid minimal worksheet (hamil + 1-det multidet)."""
    return make_worksheet(multidet=[(0, 1.0, 0.0)])


def _assert_one(run_ph2, path, fragment):
    rc, out, _ = run_ph2("validate", path)
    assert rc == 1
    line = out.strip().splitlines()[0]
    assert line.startswith(f"{path}: /")
    assert fragment in line


def test_ws_baseline_clean(run_ph2, ws):
    rc, out, _ = run_ph2("validate", ws)
    assert (rc, out) == (0, "")


def test_pauli_byte_out_of_range(run_ph2, ws):
    with _edit(ws) as f:
        f["pauli_hamil/paulis"][0, 0] = 4
    _assert_one(run_ph2, ws, "byte outside {0,1,2,3}")


def test_identity_row(run_ph2, ws):
    with _edit(ws) as f:
        f["pauli_hamil/paulis"][0, :] = 0
    _assert_one(run_ph2, ws, "all-identity row")


def test_coeffs_bad_shape(run_ph2, ws):
    with _edit(ws) as f:
        del f["state_prep/multidet/coeffs"]
        f["state_prep/multidet"].create_dataset(
            "coeffs", data=np.zeros((1, 3)))
    _assert_one(run_ph2, ws, "not float64 (M, 2)")


def test_dets_byte_out_of_range(run_ph2, ws):
    with _edit(ws) as f:
        f["state_prep/multidet/dets"][0, 0] = 2
    _assert_one(run_ph2, ws, "byte outside {0,1}")


def test_missing_normalization(run_ph2, ws):
    with _edit(ws) as f:
        del f["pauli_hamil"].attrs["normalization"]
    _assert_one(run_ph2, ws, "normalization")


def test_nonpositive_normalization(run_ph2, ws):
    with _edit(ws) as f:
        f["pauli_hamil"].attrs["normalization"] = np.float64(0.0)
    _assert_one(run_ph2, ws, "normalization")


def test_missing_state_prep(run_ph2, make_worksheet):
    pak = make_worksheet()
    _assert_one(run_ph2, pak, "/state_prep: missing")


def test_hamil_only_accepts_pak(run_ph2, make_worksheet):
    pak = make_worksheet()
    rc, out, _ = run_ph2("validate", "--hamil-only", pak)
    assert (rc, out) == (0, "")


def test_multidet_nqb_mismatch(run_ph2, make_worksheet):
    ws = make_worksheet(multidet=[(0, 1.0, 0.0)])
    with _edit(ws) as f:
        del f["state_prep/multidet/dets"]
        f["state_prep/multidet"].create_dataset(
            "dets", data=np.zeros((1, 3), dtype=np.uint8))
    _assert_one(run_ph2, ws, "nqb differs")


def _coeff_ws(make_worksheet, **attrs):
    """hamil + coeff_matrix worksheet; attrs override defaults."""
    ws = make_worksheet(nqb=8)
    a = {"n_qubits": 8, "n_sites": 4, "n_alpha": 2, "n_beta": 2,
         "closed_shell": 1, "tapered": 0}
    a.update(attrs)
    with _edit(ws) as f:
        cm = f.create_group("state_prep/coeff_matrix")
        cm.attrs["n_qubits"] = np.uint32(a["n_qubits"])
        cm.attrs["n_sites"] = np.uint32(a["n_sites"])
        cm.attrs["n_alpha"] = np.uint32(a["n_alpha"])
        cm.attrs["n_beta"] = np.uint32(a["n_beta"])
        cm.attrs["closed_shell"] = np.uint8(a["closed_shell"])
        cm.attrs["tapered"] = np.uint8(a["tapered"])
        cm.create_dataset("C_alpha", data=np.eye(4)[:, :2])
    return ws


def test_coeff_baseline_clean(run_ph2, make_worksheet):
    ws = _coeff_ws(make_worksheet)
    rc, out, _ = run_ph2("validate", ws)
    assert (rc, out) == (0, "")


def test_coeff_nqubits_rule(run_ph2, make_worksheet):
    ws = _coeff_ws(make_worksheet, n_qubits=7)
    _assert_one(run_ph2, ws, "2*n_sites - 2*tapered")


def test_coeff_cbeta_on_closed(run_ph2, make_worksheet):
    ws = _coeff_ws(make_worksheet)
    with _edit(ws) as f:
        f["state_prep/coeff_matrix"].create_dataset(
            "C_beta", data=np.eye(4)[:, :2])
    _assert_one(run_ph2, ws, "present on closed shell")


def test_coeff_cbeta_missing_on_open(run_ph2, make_worksheet):
    ws = _coeff_ws(make_worksheet, closed_shell=0)
    _assert_one(run_ph2, ws, "missing on open shell")


def test_csf_subgroup_numbering(run_ph2, make_worksheet):
    ws = _coeff_ws(make_worksheet)
    with _edit(ws) as f:
        csf = f["state_prep/coeff_matrix"].create_group("csf")
        csf.attrs["n_components"] = np.uint32(2)
        for name in ("0", "2"):
            blk = csf.create_group(name)
            blk.attrs["coefficient"] = np.float64(0.5)
            blk.create_dataset("C_alpha", data=np.eye(4)[:, :2])
    _assert_one(run_ph2, ws, "subgroups not named")


def test_csf_missing_coefficient(run_ph2, make_worksheet):
    ws = _coeff_ws(make_worksheet)
    with _edit(ws) as f:
        csf = f["state_prep/coeff_matrix"].create_group("csf")
        csf.attrs["n_components"] = np.uint32(1)
        blk = csf.create_group("0")
        blk.create_dataset("C_alpha", data=np.eye(4)[:, :2])
    _assert_one(run_ph2, ws, "coefficient")


def test_nan_before_real_row(run_ph2, make_worksheet):
    v = np.zeros((4, 2))
    v[1] = np.nan
    ws = make_worksheet(multidet=[(0, 1.0, 0.0)],
                        results={"circ_trott": ({"delta": 0.1}, v)})
    _assert_one(run_ph2, ws, "NaN row precedes")


def test_missing_delta(run_ph2, make_worksheet):
    ws = make_worksheet(multidet=[(0, 1.0, 0.0)],
                        results={"circ_trott": ({}, np.zeros((2, 2)))})
    _assert_one(run_ph2, ws, "@delta: missing")


def test_qdrift_attrs(run_ph2, make_worksheet):
    base = {"step_size": 0.1, "num_samples": 2, "depth": 4}
    v = np.zeros((2, 2))
    # without seed: passes (optional)
    ws = make_worksheet(name="a.h5", multidet=[(0, 1.0, 0.0)],
                        results={"circ_qdrift": (base, v)})
    rc, out, _ = run_ph2("validate", ws)
    assert (rc, out) == (0, "")
    # missing num_samples: fails
    bad = {k: v_ for k, v_ in base.items() if k != "num_samples"}
    ws = make_worksheet(name="b.h5", multidet=[(0, 1.0, 0.0)],
                        results={"circ_qdrift": (bad, v)})
    _assert_one(run_ph2, ws, "@num_samples: missing")


def test_cmpsit_missing_seed(run_ph2, make_worksheet):
    attrs = {"length": 1, "depth": 1, "angle_det": 0.1,
             "angle_rand": 0.1, "steps": 1}
    ws = make_worksheet(multidet=[(0, 1.0, 0.0)],
                        results={"circ_cmpsit": (attrs,
                                                 np.zeros((1, 2)))})
    _assert_one(run_ph2, ws, "@seed: missing")


def test_multiple_files(run_ph2, testdata):
    rc, out, _ = run_ph2("validate", testdata / "N4_both.h5",
                         testdata / "N4_csf_empty.h5")
    assert rc == 1
    assert "N4_both.h5" in out and "N4_csf_empty.h5" in out
