#!/usr/bin/env python3

"""ph2 -- inspect, validate, build, analyse and edit phase2 worksheets.

A worksheet (simul.h5) carries a Pauli Hamiltonian, a state-prep
payload and per-algorithm result groups; the on-disk layout is
specified in doc/simul-h5-specs.md and documented for this tool in
doc/ph2.md.  ph2 is the user-facing toolkit around that file format:

    show        summarise worksheet contents
    validate    check worksheets against the specification

Exit codes: 0 success; 1 semantic outcome (violations found, files
differ); 2 usage, IO or missing-dependency error.  Diagnostics go to
stderr prefixed "ph2: "; machine-readable output goes to stdout.

This driver only dispatches; the implementation lives in the private
package util/_ph2.  Subcommands import their modules lazily, so
`ph2.py --help` works without h5py installed.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))


def cmd_show(args):
    from _ph2 import worksheet
    return worksheet.cmd_show(args)


def cmd_validate(args):
    from _ph2 import validate
    return validate.cmd_validate(args)


def make_parser():
    from _ph2 import __version__  # stdlib-only import
    p = argparse.ArgumentParser(
        prog="ph2",
        description="Inspect, validate, build, analyse and edit"
                    " phase2 worksheets (simul.h5).",
        epilog="Full documentation: doc/ph2.md",
    )
    p.add_argument("--version", action="version",
                   version=f"ph2 {__version__}")
    sub = p.add_subparsers(dest="command", required=True,
                           metavar="COMMAND")

    q = sub.add_parser("show", help="summarise worksheet contents")
    q.add_argument("files", nargs="+", metavar="FILE")
    q.set_defaults(func=cmd_show)

    q = sub.add_parser(
        "validate",
        help="check worksheets against doc/simul-h5-specs.md")
    q.add_argument("files", nargs="+", metavar="FILE")
    q.add_argument("--hamil-only", action="store_true",
                   help="do not require /state_prep (intermediate"
                        " paks from `ph2 hamil`)")
    q.set_defaults(func=cmd_validate)

    return p


def main(argv=None):
    args = make_parser().parse_args(argv)
    from _ph2 import Ph2Error, Ph2Failure
    try:
        return args.func(args)
    except Ph2Failure as e:
        print(f"ph2: {e}", file=sys.stderr)
        return 1
    except Ph2Error as e:
        print(f"ph2: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
