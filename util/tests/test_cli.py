"""Driver-level conventions and ph2 show (CLI specification)."""

import numpy as np
import pytest

pytest.importorskip("h5py")


# -- global conventions ------------------------------------------------

def test_help_exits_zero(run_ph2):
    rc, out, _ = run_ph2("--help")
    assert rc == 0
    assert "COMMAND" in out


def test_subcommand_help_exits_zero(run_ph2):
    rc, _, _ = run_ph2("show", "--help")
    assert rc == 0


def test_version(run_ph2):
    rc, out, _ = run_ph2("--version")
    assert rc == 0
    assert out.startswith("ph2 ")


def test_unknown_subcommand_exit2(run_ph2):
    rc, out, err = run_ph2("frobnicate")
    assert rc == 2
    assert out == ""
    assert err != ""


def test_missing_file_exit2(run_ph2, tmp_path):
    rc, _, err = run_ph2("show", tmp_path / "absent.h5")
    assert rc == 2
    assert "ph2: " in err


# -- ph2 show ----------------------------------------------------------

def test_show_coeff_matrix_fixture(run_ph2, testdata):
    rc, out, _ = run_ph2("show", testdata / "N4_closed.h5")
    assert rc == 0
    assert ("coeff_matrix: n_sites=4 n_alpha=2 n_beta=2"
            " closed tapered=0") in out


def test_show_csf_components(run_ph2, testdata):
    rc, out, _ = run_ph2("show", testdata / "N4_csf.h5")
    assert rc == 0
    assert "csf: 2 components" in out


def test_show_hamil_and_multidet(run_ph2, examples_data):
    rc, out, _ = run_ph2("show", examples_data / "hamil.h5")
    assert rc == 0
    assert "251 terms, 10 qubits" in out
    assert "multidet: 1 det, 10 qubits" in out


def test_show_progress_counter(run_ph2, make_worksheet):
    values = np.full((8, 2), np.nan)
    values[:3] = 0.5
    ws = make_worksheet(
        results={"circ_trott": ({"delta": 0.1}, values)})
    rc, out, _ = run_ph2("show", ws)
    assert rc == 0
    assert "delta=0.1, 3/8 steps" in out


def test_show_qdrift_cmpsit_attrs(run_ph2, make_worksheet):
    v = np.zeros((2, 2))
    ws = make_worksheet(results={
        "circ_qdrift": ({"step_size": 0.125, "depth": 8,
                         "num_samples": 2, "seed": 7}, v),
        "circ_cmpsit": ({"length": 4, "depth": 2,
                         "angle_det": 0.0625, "angle_rand": 0.0625,
                         "steps": 8, "seed": 9}, v),
    })
    rc, out, _ = run_ph2("show", ws)
    assert rc == 0
    assert "step_size=0.125 depth=8 seed=7, 2/2 samples" in out
    assert ("length=4 depth=2 angle_det=0.0625 angle_rand=0.0625"
            " steps=8 seed=9, 2/2 samples") in out


def test_show_unrecognised_group(run_ph2, make_worksheet):
    import h5py
    ws = make_worksheet()
    with h5py.File(ws, "a") as f:
        f.create_group("foo")
    rc, out, _ = run_ph2("show", ws)
    assert rc == 0
    assert "/foo" in out
    assert "(unrecognised)" in out


def test_show_multiple_files_continue(run_ph2, make_worksheet,
                                      tmp_path):
    ws = make_worksheet()
    rc, out, err = run_ph2("show", tmp_path / "absent.h5", ws)
    assert rc == 2
    assert f"== {ws} ==" in out
    assert "pauli_hamil" in out
    assert "ph2: " in err
