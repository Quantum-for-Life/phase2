"""ph2 energy: fft, mc, ref (rpe in its own section below)."""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

TAU = 2 * np.pi


def _csv3(out):
    e, e_ref, de = (float(x) for x in out.strip().split(","))
    return e, e_ref, de


def _tone(n, phase_per_step):
    """values[s] = exp(i * phase * (s+1)) as an (n, 2) array."""
    s = np.arange(1, n + 1)
    z = np.exp(1j * phase_per_step * s)
    return np.column_stack([z.real, z.imag])


# -- output contract / E_ref --------------------------------------------

def test_ref_regression(run_ph2, examples_data):
    rc, out, _ = run_ph2("energy", "ref", examples_data / "hamil.h5")
    assert rc == 0
    assert float(out) == pytest.approx(-74.963147, abs=2e-6)


def test_ref_multidet_equals_coeff(run_ph2, testdata):
    # Same physical state, same tiny fixture Hamiltonian: the two
    # state-prep code paths must agree.
    _, out1, _ = run_ph2("energy", "ref", testdata / "N4_closed.h5")
    _, out2, _ = run_ph2("energy", "ref",
                         testdata / "N4_multidet.h5")
    assert float(out1) == pytest.approx(float(out2), abs=1e-9)


def test_eref_nan_without_state_prep(run_ph2, make_worksheet):
    ws = make_worksheet(
        results={"circ_trott": ({"delta": 0.1}, _tone(8, 0.05))})
    rc, out, _ = run_ph2("energy", "fft", ws)
    assert rc == 0
    e, e_ref, de = out.strip().split(",")
    assert e_ref == "nan" and de == "nan"


def test_missing_group_exit2(run_ph2, make_worksheet):
    ws = make_worksheet(multidet=[(0, 1.0, 0.0)])
    rc, _, err = run_ph2("energy", "fft", ws)
    assert rc == 2
    assert "circ_trott" in err


# -- fft -----------------------------------------------------------------

def test_fft_on_bin_tone(run_ph2, make_worksheet):
    n, delta, norm, offset = 64, 0.2, 0.5, -1.0
    f_bin = -5 / n                       # negative-energy line
    lam = TAU * f_bin / delta            # normalised eigenvalue
    ws = make_worksheet(
        multidet=[(0, 1.0, 0.0)],
        results={"circ_trott": ({"delta": delta},
                                _tone(n, lam * delta))})
    with h5py.File(ws, "a") as f:
        f["pauli_hamil"].attrs["normalization"] = np.float64(norm)
        f["pauli_hamil"].attrs["offset"] = np.float64(offset)
    rc, out, _ = run_ph2("energy", "fft", ws)
    assert rc == 0
    e, _, _ = _csv3(out)
    assert e == pytest.approx(lam / norm + offset, abs=1e-9)


def test_fft_dc_excluded(run_ph2, make_worksheet):
    # A large constant bias must not win over the physical tone.
    n, delta = 64, 0.2
    lam = TAU * (-5 / n) / delta
    values = _tone(n, lam * delta)
    values[:, 0] += 3.0
    ws = make_worksheet(
        multidet=[(0, 1.0, 0.0)],
        results={"circ_trott": ({"delta": delta}, values)})
    rc, out, _ = run_ph2("energy", "fft", ws)
    norm = 1.0 / 0.75  # make_worksheet default
    e, _, _ = _csv3(out)
    assert e == pytest.approx(lam / norm + (-1.0), abs=1e-9)


def test_fft_peaks_lists_tone(run_ph2, make_worksheet):
    pytest.importorskip("scipy")
    n, delta = 64, 0.2
    lam = TAU * (-5 / n) / delta
    ws = make_worksheet(
        multidet=[(0, 1.0, 0.0)],
        results={"circ_trott": ({"delta": delta},
                                _tone(n, lam * delta))})
    rc, out, _ = run_ph2("energy", "fft", "--peaks", ws)
    assert rc == 0
    lines = out.strip().splitlines()
    assert len(lines) >= 2
    dominant = _csv3(lines[0])[0]
    peak_energies = [float(l.split(",")[0]) for l in lines[1:]]
    assert any(abs(p - dominant) < 1e-9 for p in peak_energies)


def test_fft_peaks_missing_scipy_hint(run_ph2, make_worksheet,
                                      tmp_path):
    import os
    ws = make_worksheet(
        multidet=[(0, 1.0, 0.0)],
        results={"circ_trott": ({"delta": 0.2}, _tone(8, 0.1))})
    block = tmp_path / "block"
    block.mkdir()
    (block / "scipy.py").write_text(
        'raise ImportError("blocked for test")\n')
    env = dict(os.environ, PYTHONPATH=str(block))
    rc, _, err = run_ph2("energy", "fft", "--peaks", ws, env=env)
    assert rc == 2
    assert "[examples]" in err


# -- mc ------------------------------------------------------------------

def _mc_ws(make_worksheet, group, attrs, phase, n=16, pad_nan=0):
    z = np.exp(1j * phase)
    values = np.tile([z.real, z.imag], (n, 1))
    if pad_nan:
        values = np.vstack([values, np.full((pad_nan, 2), np.nan)])
    return make_worksheet(multidet=[(0, 1.0, 0.0)],
                          results={group: (attrs, values)})


def test_mc_qdrift_formula(run_ph2, make_worksheet):
    attrs = {"step_size": 0.125, "depth": 8, "num_samples": 16,
             "seed": 1}
    phase = -0.3
    ws = _mc_ws(make_worksheet, "circ_qdrift", attrs, phase)
    rc, out, _ = run_ph2("energy", "mc", ws)
    assert rc == 0
    T = 8 * np.arcsin(0.125)
    norm, offset = 1.0 / 0.75, -1.0
    e, _, _ = _csv3(out)
    assert e == pytest.approx(phase / (T * norm) + offset,
                              abs=1e-12)


def test_mc_cmpsit_formula_and_nan_suffix(run_ph2, make_worksheet):
    attrs = {"length": 4, "depth": 2, "angle_det": 0.0625,
             "angle_rand": 0.0625, "steps": 8, "seed": 1}
    phase = -0.2
    ws = _mc_ws(make_worksheet, "circ_cmpsit", attrs, phase,
                pad_nan=4)
    rc, out, _ = run_ph2("energy", "mc", "--group", "circ_cmpsit",
                         ws)
    assert rc == 0
    T = 8 * 0.0625
    norm, offset = 1.0 / 0.75, -1.0
    e, _, _ = _csv3(out)
    assert e == pytest.approx(phase / (T * norm) + offset,
                              abs=1e-12)


def test_mc_decoherence_warning(run_ph2, make_worksheet):
    values = np.array([[1.0, 0.0], [-1.0, 0.0]] * 8, dtype=float)
    attrs = {"step_size": 0.125, "depth": 8, "num_samples": 16,
             "seed": 1}
    ws = make_worksheet(multidet=[(0, 1.0, 0.0)],
                        results={"circ_qdrift": (attrs, values)})
    rc, _, err = run_ph2("energy", "mc", ws)
    assert rc == 0
    assert "decohered" in err
