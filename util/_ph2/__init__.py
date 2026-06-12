"""Private internals of util/ph2.py.

Not a user-accessible library: the only supported surface is the
ph2.py command line.  This package must import without h5py/numpy so
that `ph2.py --help` works in a bare interpreter; submodules import
the heavy dependencies at their own top level instead.
"""

__version__ = "1.1.0"

# Result groups a worksheet may carry, with the attributes each one
# is required to have (doc/simul-h5-specs.md; /circ_qdrift "seed" is
# written by current ph2run but optional for older files).
RESULTS_GROUPS = {
    "circ_trott": ("delta",),
    "circ_trott2": ("delta",),
    "circ_qdrift": ("step_size", "num_samples", "depth"),
    "circ_cmpsit": ("length", "depth", "angle_det", "angle_rand",
                    "steps", "seed"),
}


class Ph2Error(Exception):
    """Usage / IO / missing-dependency error: exit code 2."""


class Ph2Failure(Exception):
    """Semantic outcome reported as failure: exit code 1."""


def open_ro(path):
    """Open an HDF5 file read-only; Ph2Error on any failure."""
    import h5py
    try:
        return h5py.File(path, "r")
    except OSError as e:
        raise Ph2Error(f"cannot open {path}: {e}") from None


def guard_out(path, force):
    """Refuse to overwrite an existing output unless --force."""
    import os
    if os.path.exists(path) and not force:
        raise Ph2Error(f"{path} exists (use --force to overwrite)")
