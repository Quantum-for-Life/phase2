#!/usr/bin/env python

import argparse
import h5py
import shutil


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="qdrift_sample",
        description="Sample Hamiltonian for QDRIFT",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("filename", type=str, help="Input file")
    parser.add_argument("--num-samples", type=int, default=128, help="number of QDIRFT samples")
    parser.add_argument("--step-size", type=float)
    parser.add_argument("--depth", type=int)
    parser.add_argument("--outs", type=int)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    for i in range(1, args.outs + 1):
        filename = args.filename + f"-{i:04}"
        shutil.copy(args.filename, filename)
        with h5py.File(filename, "a") as f:
            grp = f.create_group("trotter_steps")
            grp.attrs["step_size"] = args.step_size
            grp.attrs["num_samples"] = args.num_samples
            grp.attrs["depth"] = args.depth
