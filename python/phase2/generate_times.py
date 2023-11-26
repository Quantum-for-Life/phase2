#!/bin/env/python

import argparse

import h5py
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="jw_map",
        description="Parse FCIDUMP to JW mapping",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("-o", "--output", type=str, help="Output file: .h5")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


STEPS = 111
TIMES = [float(x) for x in range(0, STEPS)]
VALUES = np.ndarray((STEPS, 2), dtype='d')

for i in range(0, STEPS):
    VALUES[i, 0] = np.nan
    VALUES[i, 1] = np.nan


def h5_output(outfile: str):
    with h5py.File(outfile, "a") as f:
        grp = f.create_group("time_series")
        grp.create_dataset("times", (STEPS,), dtype='d')[...] = TIMES
        grp.create_dataset("values", (STEPS, 2), dtype='d')[...] = VALUES


if __name__ == "__main__":
    args = parse_arguments()

    if args.output:
        h5_output(args.output)
