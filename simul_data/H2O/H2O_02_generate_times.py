#!/bin/env/python

import argparse
import random
import shutil
from pathlib import Path

import h5py
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("-n", "--out-files", type=int, default=1,
                        help="Number of output files")
    parser.add_argument("-s", "--time-steps", type=int, default=128,
                        help="Number of time steps in each file")
    parser.add_argument("--time-min", type=float, default=0.0)
    parser.add_argument("--time-max", type=float, default=100.0)
    return parser.parse_args()


def write_h5file(filename, times):
    size = len(times)
    values = np.ndarray((size, 2), dtype="d")
    for i in range(0, size):
        values[i, 0] = np.nan
        values[i, 1] = np.nan
    with h5py.File(filename, "a") as f:
        grp = f.create_group("time_series")
        grp.create_dataset("times", (size,), dtype="d")[...] = times
        grp.create_dataset("values", (size, 2), dtype="d")[...] = values


if __name__ == "__main__":
    args = parse_arguments()

    num_files = args.out_files
    time_steps = args.time_steps
    time_min = args.time_min
    time_max = args.time_max

    times = set()
    while len(times) < num_files * time_steps:
        times.add(random.randint(time_min, time_max))
    times = list(times)
    random.shuffle(times)

    h5file = args.filename
    for fi in range(args.out_files):
        filename = Path(h5file).stem + f"_{fi:03}.h5";
        shutil.copy(h5file, filename)
        time_slice = times[fi*time_steps:(fi+1)*time_steps]
        write_h5file(filename, time_slice)
