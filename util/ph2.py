#!/usr/bin/env python3

"""ph2 -- inspect, validate, build, analyse and edit phase2 worksheets.

A worksheet (simul.h5) carries a Pauli Hamiltonian, a state-prep
payload and per-algorithm result groups; the on-disk layout is
specified in doc/simul-h5-specs.md and documented for this tool in
doc/ph2.md.  ph2 is the user-facing toolkit around that file format:

    show        summarise worksheet contents
    validate    check worksheets against the specification
    hamil       build /pauli_hamil (FCIDUMP or Pauli-term text)
    stprep      build /state_prep (multidet or coeff_matrix)
    energy      extract energies (fft, mc, ref, rpe, rpe-qdrift)

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


def cmd_hamil_paulis(args):
    from _ph2 import hamil
    return hamil.cmd_paulis(args)


def cmd_hamil_fcidump(args):
    from _ph2 import hamil
    return hamil.cmd_fcidump(args)


def cmd_stprep_multidet(args):
    from _ph2 import stprep
    return stprep.cmd_multidet(args)


def cmd_stprep_coeff(args):
    from _ph2 import stprep
    return stprep.cmd_coeff(args)


def cmd_energy_fft(args):
    from _ph2 import analysis
    return analysis.cmd_fft(args)


def cmd_energy_mc(args):
    from _ph2 import analysis
    return analysis.cmd_mc(args)


def cmd_energy_ref(args):
    from _ph2 import analysis
    return analysis.cmd_ref(args)


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

    q = sub.add_parser("hamil", help="build /pauli_hamil")
    hs = q.add_subparsers(dest="builder", required=True,
                          metavar="BUILDER")
    b = hs.add_parser("paulis",
                      help="from a Pauli-term text file (TERMS)")
    b.add_argument("terms", metavar="TERMS")
    b.add_argument("-o", dest="out", required=True,
                   metavar="OUT.h5")
    b.add_argument("--n-qubits", type=int, metavar="N",
                   help="widen the register beyond max index + 1")
    b.add_argument("--offset", type=float, default=0.0,
                   metavar="E0",
                   help="extra identity contribution (Hartree)")
    b.add_argument("--sort-terms", action="store_true",
                   help="sort terms by Pauli label")
    b.add_argument("--force", action="store_true",
                   help="overwrite an existing output file")
    b.set_defaults(func=cmd_hamil_paulis)

    b = hs.add_parser(
        "fcidump",
        help="Jordan-Wigner transform of an FCIDUMP file")
    b.add_argument("fcidump", metavar="FCIDUMP")
    b.add_argument("-o", dest="out", required=True,
                   metavar="OUT.h5")
    b.add_argument("--sort-terms", action="store_true",
                   help="sort terms by Pauli label")
    b.add_argument("--force", action="store_true",
                   help="overwrite an existing output file")
    b.set_defaults(func=cmd_hamil_fcidump)

    q = sub.add_parser("stprep", help="build /state_prep")
    ss = q.add_subparsers(dest="builder", required=True,
                          metavar="BUILDER")
    b = ss.add_parser("multidet",
                      help="from an INPUTST determinant list")
    b.add_argument("inputst", metavar="INPUTST")
    b.add_argument("-f", dest="file", required=True, metavar="FILE",
                   help="worksheet to append to")
    b.add_argument("--no-reorder", action="store_true",
                   help="occupations are already alpha-block |"
                        " beta-block")
    b.add_argument("--force", action="store_true",
                   help="replace an existing /state_prep subtype")
    b.set_defaults(func=cmd_stprep_multidet)

    b = ss.add_parser("coeff",
                      help="from CMAT coefficient-matrix files")
    b.add_argument("-f", dest="file", required=True, metavar="FILE",
                   help="worksheet to append to")
    b.add_argument("--c-alpha", required=True, metavar="CA",
                   help="alpha-block CMAT file")
    b.add_argument("--c-beta", metavar="CB",
                   help="beta-block CMAT file (open shell)")
    b.add_argument("--tapered", action="store_true",
                   help="Z2 + S_z tapering"
                        " (n_qubits = 2*n_sites - 2)")
    b.add_argument("--csf", action="append", metavar="W:CA[,CB]",
                   help="CSF component (repeatable, ordered)")
    b.add_argument("--force", action="store_true",
                   help="replace an existing /state_prep subtype")
    b.set_defaults(func=cmd_stprep_coeff)

    q = sub.add_parser("energy",
                       help="extract energies from results")
    es = q.add_subparsers(dest="method", required=True,
                          metavar="METHOD")
    b = es.add_parser("fft",
                      help="FFT of a Trotter overlap time series")
    b.add_argument("filename", metavar="FILE")
    b.add_argument("--group", default="circ_trott",
                   choices=("circ_trott", "circ_trott2"))
    b.add_argument("--peaks", action="store_true",
                   help="also list every detected spectral peak")
    b.set_defaults(func=cmd_energy_fft)

    b = es.add_parser("mc",
                      help="phase of the averaged Monte-Carlo"
                           " overlap")
    b.add_argument("filename", metavar="FILE")
    b.add_argument("--group", default="circ_qdrift",
                   choices=("circ_qdrift", "circ_cmpsit"))
    b.set_defaults(func=cmd_energy_mc)

    b = es.add_parser("ref",
                      help="trial-state energy <psi|H|psi>")
    b.add_argument("filename", metavar="FILE")
    b.set_defaults(func=cmd_energy_ref)

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
