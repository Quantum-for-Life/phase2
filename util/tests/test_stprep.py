"""ph2 stprep: multidet and coeff builders, expansion oracle."""

import importlib.util

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")


def _savetxt(path, mat):
    np.savetxt(path, mat)
    return path


def _hamil_pak(run_ph2, tmp_path, nqb):
    """Hamil-only pak of the requested width via ph2 itself."""
    terms = tmp_path / f"T{nqb}"
    terms.write_text(f"1.0 Z0 Z{nqb - 1}\n")
    pak = tmp_path / f"pak{nqb}.h5"
    rc, _, err = run_ph2("hamil", "paulis", terms, "-o", pak)
    assert rc == 0, err
    return pak


# -- multidet ----------------------------------------------------------

def test_multidet_round_trip(run_ph2, examples_data, tmp_path):
    pak = _hamil_pak(run_ph2, tmp_path, 10)
    rc, _, err = run_ph2("stprep", "multidet",
                         examples_data / "INPUTST", "-f", pak)
    assert rc == 0, err
    with h5py.File(pak) as f, \
            h5py.File(examples_data / "hamil.h5") as ref:
        assert np.array_equal(f["state_prep/multidet/coeffs"][...],
                              ref["state_prep/multidet/coeffs"][...])
        assert np.array_equal(f["state_prep/multidet/dets"][...],
                              ref["state_prep/multidet/dets"][...])
    rc, out, _ = run_ph2("validate", pak)
    assert (rc, out) == (0, "")


def test_multidet_no_reorder(run_ph2, examples_data, tmp_path):
    a = _hamil_pak(run_ph2, tmp_path, 10)
    run_ph2("stprep", "multidet", examples_data / "INPUTST",
            "-f", a, "--no-reorder")
    with h5py.File(a) as f, \
            h5py.File(examples_data / "hamil.h5") as ref:
        assert not np.array_equal(
            f["state_prep/multidet/dets"][...],
            ref["state_prep/multidet/dets"][...])


def test_multidet_nqb_mismatch_exit1(run_ph2, examples_data,
                                     tmp_path):
    pak = _hamil_pak(run_ph2, tmp_path, 4)
    rc, _, err = run_ph2("stprep", "multidet",
                         examples_data / "INPUTST", "-f", pak)
    assert rc == 1
    assert "qubits" in err


def test_stprep_ambiguity_guard(run_ph2, examples_data, tmp_path):
    pak = _hamil_pak(run_ph2, tmp_path, 10)
    run_ph2("stprep", "multidet", examples_data / "INPUTST",
            "-f", pak)
    rc, _, err = run_ph2("stprep", "multidet",
                         examples_data / "INPUTST", "-f", pak)
    assert rc == 1
    assert "already carries" in err
    rc, _, err = run_ph2("stprep", "multidet",
                         examples_data / "INPUTST", "-f", pak,
                         "--force")
    assert rc == 0, err
    rc, out, _ = run_ph2("validate", pak)
    assert (rc, out) == (0, "")


# -- coeff: fixture round-trips ----------------------------------------

def _compare_group(got, want):
    """Attrs (incl. dtypes) and datasets of two h5 groups equal."""
    assert set(got.attrs) == set(want.attrs)
    for k in want.attrs:
        assert got.attrs[k] == want.attrs[k], k
        assert np.asarray(got.attrs[k]).dtype == \
            np.asarray(want.attrs[k]).dtype, k
    assert set(got.keys()) == set(want.keys())
    for k in want:
        if isinstance(want[k], h5py.Group):
            _compare_group(got[k], want[k])
        else:
            assert np.array_equal(got[k][...], want[k][...]), k
            assert got[k].dtype == want[k].dtype, k


def _roundtrip(run_ph2, tmp_path, testdata, fixture, *flags):
    """Rebuild `fixture`'s coeff_matrix via ph2 stprep coeff."""
    with h5py.File(testdata / fixture) as f:
        cm = f["state_prep/coeff_matrix"]
        nqb = int(cm.attrs["n_qubits"])
        ca = _savetxt(tmp_path / "CA", cm["C_alpha"][...])
        args = ["--c-alpha", ca]
        if "C_beta" in cm:
            args += ["--c-beta",
                     _savetxt(tmp_path / "CB", cm["C_beta"][...])]
        if "csf" in cm:
            for k in range(int(cm["csf"].attrs["n_components"])):
                blk = cm[f"csf/{k}"]
                p = _savetxt(tmp_path / f"CSF{k}",
                             blk["C_alpha"][...])
                args += ["--csf",
                         f"{blk.attrs['coefficient']}:{p}"]
    pak = _hamil_pak(run_ph2, tmp_path, nqb)
    rc, _, err = run_ph2("stprep", "coeff", "-f", pak,
                         *args, *flags)
    assert rc == 0, err
    with h5py.File(pak) as got, h5py.File(testdata / fixture) as ref:
        _compare_group(got["state_prep/coeff_matrix"],
                       ref["state_prep/coeff_matrix"])
    rc, out, _ = run_ph2("validate", pak)
    assert (rc, out) == (0, "")


def test_coeff_closed(run_ph2, tmp_path, testdata):
    _roundtrip(run_ph2, tmp_path, testdata, "N4_closed.h5")


def test_coeff_open(run_ph2, tmp_path, testdata):
    _roundtrip(run_ph2, tmp_path, testdata, "N4_open.h5")


def test_coeff_csf(run_ph2, tmp_path, testdata):
    _roundtrip(run_ph2, tmp_path, testdata, "N4_csf.h5")


def test_coeff_tapered(run_ph2, tmp_path, testdata):
    _roundtrip(run_ph2, tmp_path, testdata, "N8_tapered.h5",
               "--tapered")


# -- CMAT grammar -------------------------------------------------------

def test_cmat_comments_ok(run_ph2, tmp_path):
    pak = _hamil_pak(run_ph2, tmp_path, 4)
    ca = tmp_path / "CA"
    ca.write_text("# comment\n 1.0 0.0\n\n 0.0 1.0  # x\n")
    rc, _, err = run_ph2("stprep", "coeff", "-f", pak,
                         "--c-alpha", ca)
    assert rc == 0, err
    with h5py.File(pak) as f:
        cm = f["state_prep/coeff_matrix"]
        assert cm.attrs["n_sites"] == 2
        assert cm.attrs["n_alpha"] == 2


def test_cmat_ragged_exit2(run_ph2, tmp_path):
    pak = _hamil_pak(run_ph2, tmp_path, 4)
    ca = tmp_path / "CA"
    ca.write_text("1.0 0.0\n0.0\n")
    rc, _, err = run_ph2("stprep", "coeff", "-f", pak,
                         "--c-alpha", ca)
    assert rc == 2
    assert "ragged" in err


def test_cmat_beta_rows_mismatch_exit2(run_ph2, tmp_path):
    pak = _hamil_pak(run_ph2, tmp_path, 4)
    ca = _savetxt(tmp_path / "CA", np.eye(2))
    cb = _savetxt(tmp_path / "CB", np.eye(3))
    rc, _, err = run_ph2("stprep", "coeff", "-f", pak,
                         "--c-alpha", ca, "--c-beta", cb)
    assert rc == 2
    assert "row count" in err


def test_csf_shape_mismatch_exit2(run_ph2, tmp_path):
    pak = _hamil_pak(run_ph2, tmp_path, 4)
    ca = _savetxt(tmp_path / "CA", np.eye(2))
    bad = _savetxt(tmp_path / "BAD", np.eye(3))
    rc, _, err = run_ph2("stprep", "coeff", "-f", pak,
                         "--c-alpha", ca, "--csf", f"0.5:{bad}")
    assert rc == 2
    assert "shape" in err


def test_csf_shell_mismatch_exit2(run_ph2, tmp_path):
    pak = _hamil_pak(run_ph2, tmp_path, 4)
    ca = _savetxt(tmp_path / "CA", np.eye(2))
    rc, _, err = run_ph2("stprep", "coeff", "-f", pak,
                         "--c-alpha", ca, "--csf",
                         f"0.5:{ca},{ca}")
    assert rc == 2
    assert "shell" in err


# -- expansion oracle ----------------------------------------------------

def _reference_oracle(repo):
    spec = importlib.util.spec_from_file_location(
        "coeff_matrix_reference",
        repo / "test" / "ref" / "coeff_matrix_reference.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.expand_coeff_matrix


@pytest.mark.parametrize("fixture,tapered",
                         [("N4_closed.h5", False),
                          ("N8_tapered.h5", True)])
def test_oracle_parity(repo, testdata, fixture, tapered):
    from _ph2.stprep import expand_coeff_matrix
    ref_expand = _reference_oracle(repo)
    with h5py.File(testdata / fixture) as f:
        cm = f["state_prep/coeff_matrix"]
        ns, na, nb = (int(cm.attrs[k])
                      for k in ("n_sites", "n_alpha", "n_beta"))
        C = cm["C_alpha"][...]
    got = expand_coeff_matrix(ns, na, nb, C, None, tapered)
    want = ref_expand(ns, na, nb, C, None, tapered)
    assert set(got) == set(want)
    for k in want:
        assert got[k] == pytest.approx(want[k], abs=1e-12)


def test_expand_matches_multidet_fixture(testdata):
    from _ph2.stprep import expand_coeff_matrix
    with h5py.File(testdata / "N4_closed.h5") as f:
        cm = f["state_prep/coeff_matrix"]
        amps = expand_coeff_matrix(
            int(cm.attrs["n_sites"]), int(cm.attrs["n_alpha"]),
            int(cm.attrs["n_beta"]), cm["C_alpha"][...], None,
            False)
    with h5py.File(testdata / "N4_multidet.h5") as f:
        coeffs = f["state_prep/multidet/coeffs"][...]
        dets = f["state_prep/multidet/dets"][...]
    want = {}
    for (re_, _), occ in zip(coeffs, dets):
        idx = sum(int(b) << j for j, b in enumerate(occ))
        want[idx] = re_
    assert set(amps) == set(want)
    for k in want:
        assert amps[k] == pytest.approx(want[k], abs=1e-12)
