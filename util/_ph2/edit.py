"""Worksheet editing (ph2 strip / diff / attr).

strip removes results groups by copy-rewrite: h5py's `del` only
unlinks and HDF5 never returns the freed pages without a repack, so
the result is written to a fresh file (-o OUT) or to a same-
directory temp file that atomically replaces the original.  Root
attributes survive; dangling legacy links are dropped (they cannot
be copied).  The in-place path breaks hard links by construction.

diff walks groups, datasets and attributes of two worksheets and
prints one difference per line.  Floats compare within an absolute
tolerance with NaN == NaN (so identically NaN-padded values
compare clean); integers and strings compare exactly; the root
`uuid` is ignored unless --strict, since two builds of the same
input differ only there.

attr reads or writes scalar attributes with an explicit dtype, so
the spec's HDF5 types (double, unsigned long, the uint8 flags) are
preservable from the command line.  This is the sharp tool: run
`ph2 validate` afterwards.
"""

import sys

import numpy as np

from . import (Ph2Error, Ph2Failure, RESULTS_GROUPS, copy_without,
               guard_out, open_ro)


def cmd_strip(args):
    with open_ro(args.file) as f:
        present = [n for n in RESULTS_GROUPS if n in f]
    if args.group:
        for name in args.group:
            if name not in present:
                raise Ph2Failure(f"{args.file}: no /{name} group")
        drop = set(args.group)
    else:
        drop = set(present)
        if not drop:
            return 0  # nothing to strip
    if args.out:
        guard_out(args.out, args.force)
        copy_without(args.file, drop, out=args.out)
    else:
        copy_without(args.file, drop)
    return 0


def _collect(f):
    """(objects, attrs): path -> h5 object / (path, name) -> value."""
    import h5py
    objects = {"/": f}
    attrs = {("/", k): v for k, v in f.attrs.items()}

    def visit(name, obj):
        path = "/" + name
        objects[path] = obj
        for k, v in obj.attrs.items():
            attrs[(path, k)] = v

    f.visititems(visit)
    return objects, attrs


def _is_results(path):
    top = path.split("/", 2)[1] if path != "/" else ""
    return top in RESULTS_GROUPS


def _values_equal(a, b, tol):
    a, b = np.asarray(a), np.asarray(b)
    if a.dtype.kind != b.dtype.kind:
        return False
    if a.shape != b.shape:
        return False
    if a.dtype.kind == "f":
        nan_a, nan_b = np.isnan(a), np.isnan(b)
        if not np.array_equal(nan_a, nan_b):
            return False
        return bool(np.all(np.abs(a[~nan_a] - b[~nan_b]) <= tol))
    return np.array_equal(a, b)


def cmd_diff(args):
    import h5py
    diffs = []
    with open_ro(args.a) as fa, open_ro(args.b) as fb:
        obj_a, attr_a = _collect(fa)
        obj_b, attr_b = _collect(fb)

        def skip(path):
            return args.ignore_results and _is_results(path)

        for path in sorted(set(obj_a) | set(obj_b)):
            if skip(path):
                continue
            if path not in obj_b:
                diffs.append(f"only in A: {path}")
                continue
            if path not in obj_a:
                diffs.append(f"only in B: {path}")
                continue
            a, b = obj_a[path], obj_b[path]
            if isinstance(a, h5py.Group) != isinstance(b, h5py.Group):
                diffs.append(f"object kind differs: {path}")
                continue
            if isinstance(a, h5py.Group):
                continue
            if a.dtype.kind != b.dtype.kind:
                diffs.append(f"dtype differs: {path}:"
                             f" {a.dtype} != {b.dtype}")
            elif a.shape != b.shape:
                diffs.append(f"shape differs: {path}:"
                             f" {a.shape} != {b.shape}")
            elif not _values_equal(a[...], b[...], args.tol):
                if a.dtype.kind == "f":
                    av, bv = a[...], b[...]
                    mask = ~np.isnan(av)
                    d = float(np.max(np.abs(av[mask] - bv[mask]),
                                     initial=0.0))
                    diffs.append(f"dataset differs: {path}:"
                                 f" max|d|={d:.3g}")
                else:
                    diffs.append(f"dataset differs: {path}")

        for key in sorted(set(attr_a) | set(attr_b)):
            path, name = key
            if skip(path):
                continue
            if not args.strict and key == ("/", "uuid"):
                continue
            ref = f"{path}@{name}"
            if key not in attr_b:
                diffs.append(f"only in A: {ref}")
                continue
            if key not in attr_a:
                diffs.append(f"only in B: {ref}")
                continue
            va, vb = attr_a[key], attr_b[key]
            ka = np.asarray(va).dtype.kind
            if ka == "f" and np.asarray(vb).dtype.kind == "f":
                same = abs(float(va) - float(vb)) <= args.tol
            else:
                same = np.array_equal(np.asarray(va),
                                      np.asarray(vb))
            if not same:
                diffs.append(f"attr differs: {ref}: {va} != {vb}")

    for line in diffs:
        print(line)
    return 1 if diffs else 0


_TYPES = {
    "f64": np.float64,
    "u64": np.uint64,
    "u32": np.uint32,
    "u8": np.uint8,
    "str": str,
}


def _resolve(f, h5path):
    try:
        return f[h5path]
    except KeyError:
        raise Ph2Error(f"no object at {h5path}") from None


def _fmt(v):
    if isinstance(v, bytes):
        return v.decode()
    return str(v)


def cmd_attr_get(args):
    with open_ro(args.file) as f:
        obj = _resolve(f, args.h5path)
        if args.name is not None:
            if args.name not in obj.attrs:
                raise Ph2Error(f"{args.h5path} has no attribute"
                               f" {args.name!r}")
            print(f"{args.name}={_fmt(obj.attrs[args.name])}")
            return 0
        for k in obj.attrs:
            print(f"{k}={_fmt(obj.attrs[k])}")
    return 0


def cmd_attr_set(args):
    import h5py
    ctor = _TYPES[args.type]
    try:
        value = ctor(args.value)
    except ValueError:
        raise Ph2Error(f"cannot parse {args.value!r} as"
                       f" {args.type}") from None
    try:
        f = h5py.File(args.file, "a")
    except OSError as e:
        raise Ph2Error(f"cannot open {args.file}: {e}") from None
    with f:
        obj = _resolve(f, args.h5path)
        if args.name in obj.attrs:
            del obj.attrs[args.name]  # allow a dtype change
        obj.attrs[args.name] = value
    return 0
