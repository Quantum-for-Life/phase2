#!/usr/bin/env python

import argparse

import h5py


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="parse_inputst",
        description="Parse INPUTST to JW mapping",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("filename", type=str, help="INPUTST input file")
    parser.add_argument("-o", "--output", type=str, help="Output file: .h5")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    with open(args.filename, 'r') as f:
        input_lines = f.read().split("\n")

    dets = []
    coeffs = []
    for line in input_lines:
        if not line:
            continue
        x = float(line.split()[0])
        y = float(line.split()[1])
        det = [int(i) for i in line.split()[2:]]
        norbs = int(len(det) / 2)
        # rearrange electrons according to our convention
        det_rearr = [det[2 * i] for i in range(norbs)] + \
                    [det[2 * i + 1] for i in range(norbs)]
        dets.append(det_rearr)
        coeffs.append(x + 1j * y)

    with h5py.File(args.output, "w") as f:

        h5_md = f.create_group("state_prep/multidet")
        h5_coeffs = h5_md.create_dataset("coeffs", shape=(len(dets), 2),
                                         dtype="d")
        h5_coeffs[...] = [[z.real, z.imag] for z in coeffs]
        h5_dets = h5_md.create_dataset("dets", shape=(len(dets), len(dets[0])),
                                       dtype="u1")
        h5_dets[...] = dets
