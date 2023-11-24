#!/bin/env/python

import argparse

import h5py
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="jw_map",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("-o", "--output", type=str, help="Output file: .h5")
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def h5_output(outfile: str):
    with h5py.File(outfile, "a") as f:
        state_prep = f.create_group("state_prep")
        multidet = state_prep.create_group("multidet")
        multidet.create_dataset("coeffs", (1,), dtype='d')[...] = [1.0]
        multidet.create_dataset("dets", (1,4), dtype='u1')[...] = [1, 0, 1, 0]


if __name__ == "__main__":
    args = parse_arguments()

    if args.output:
        h5_output(args.output)
