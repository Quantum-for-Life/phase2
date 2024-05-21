#!/usr/bin/env python

import argparse
import h5py


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="parse_fcidump",
        description="Parse FCIDUMP to JW mapping",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("-o", "--output", type=str, help="Output file: .h5")
    parser.add_argument("-t", "--time-factor", type=float,
                        help="Multiply phases by a real factor")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def h5_output(outfile: str,
              time_factor: float = 1.0):
    with h5py.File(outfile, "a") as f:
        grp = f.create_group("circ_trott")
        grp.attrs["time_factor"] = time_factor


if __name__ == "__main__":
    args = parse_arguments()

    if args.output:
        h5_output(args.output, time_factor=args.time_factor)
