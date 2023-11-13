#!/bin/env/python

import argparse

import h5py


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="jw_map",
        description="Parse FCIDUMP to JW mapping",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("-o", "--output", type=str, help="Output file: .h5")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


VALUES = [float(x) for x in range(0, 100)]

STEPS = len(VALUES)


def h5_output(outfile: str):
    with h5py.File(outfile, "a") as f:
        grp = f.create_group("time_series")
        dset_times = grp.create_dataset("times", (STEPS,), dtype='d')
        dset_times[...] = VALUES


if __name__ == "__main__":
    args = parse_arguments()

    if args.output:
        h5_output(args.output)
