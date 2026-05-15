"""Python interface to the phase2 Hamiltonian simulator.

Wraps the C function phase2_run from libphase2.so via ctypes.
The library is located by checking the PHASE2_LIB environment
variable first, then searching relative to this file.

Usage::

    import phase2

    paulis = ["ZI", "IZ", "XX"]
    coeffs = [0.5, -0.3, 0.1]
    delta = 0.01
    psi = "00"

    result = phase2.run(paulis, coeffs, delta, psi)
    print(result)  # complex number
"""

import ctypes
import os
from pathlib import Path

__all__ = ["run"]


def _load_library() -> ctypes.CDLL:
    """Locate and load libphase2.so."""
    env = os.environ.get("PHASE2_LIB")
    if env:
        p = Path(env)
        if p.is_file():
            return ctypes.CDLL(str(p))

    here = Path(__file__).resolve().parent
    candidates = [
        # In-tree development build: python/phase2/ -> repo
        # root, library lives at build/libphase2.so.
        here / ".." / ".." / "build" / "libphase2.so",
        # Legacy locations (older trees placed the .so at
        # the repo root or inside an installed package).
        here / ".." / ".." / "libphase2.so",
        here / ".." / "libphase2.so",
        here / "libphase2.so",
    ]
    for c in candidates:
        c = c.resolve()
        if c.is_file():
            return ctypes.CDLL(str(c))

    raise RuntimeError(
        "libphase2.so not found; set PHASE2_LIB or place"
        " the library alongside the package"
    )


_lib = _load_library()

_lib.phase2_run.argtypes = [
    ctypes.POINTER(ctypes.c_char_p),  # pauli_strs
    ctypes.POINTER(ctypes.c_double),  # coeffs
    ctypes.c_size_t,                  # nterms
    ctypes.c_double,                  # delta
    ctypes.c_char_p,                  # psi_str
    ctypes.POINTER(ctypes.c_double),  # out_re
    ctypes.POINTER(ctypes.c_double),  # out_im
]
_lib.phase2_run.restype = ctypes.c_int


def run(
    paulis: list[str],
    coeffs: list[float],
    delta: float,
    psi: str,
) -> complex:
    """Compute <psi| prod_k exp(i*delta*coeffs[k]*P_k) |psi>.

    Args:
        paulis: Pauli strings, e.g. ["ZI", "IZ", "XX"].
        coeffs: Real coefficients, one per Pauli string.
        delta:  Step size parameter.
        psi:    Computational-basis state as a bit string,
                e.g. "00".  The number of qubits is len(psi).

    Returns:
        Complex expectation value.

    Raises:
        ValueError: If len(paulis) != len(coeffs) or the
                    C function returns an error.
    """
    nterms = len(paulis)
    if nterms != len(coeffs):
        raise ValueError(
            f"len(paulis)={nterms} != len(coeffs)="
            f"{len(coeffs)}"
        )

    c_paulis = (ctypes.c_char_p * nterms)(
        *(s.encode("ascii") for s in paulis)
    )
    c_coeffs = (ctypes.c_double * nterms)(*coeffs)
    c_psi = psi.encode("ascii")

    out_re = ctypes.c_double(0.0)
    out_im = ctypes.c_double(0.0)

    rc = _lib.phase2_run(
        c_paulis,
        c_coeffs,
        ctypes.c_size_t(nterms),
        ctypes.c_double(delta),
        c_psi,
        ctypes.byref(out_re),
        ctypes.byref(out_im),
    )

    if rc == -1:
        raise ValueError(
            "phase2_run returned -1: invalid input"
        )

    return complex(out_re.value, out_im.value)
