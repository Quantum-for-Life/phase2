import argparse
import math


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="qdrift_sample",
        description="Compute QDrift depth and step_size from delta and epsilon",
        epilog="Quantum-for-Life",
    )

    parser.add_argument("--delta", type=float)
    parser.add_argument("--epsilon", type=float)

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    delta = args.delta
    epsilon = args.epsilon
    J = math.ceil(math.log2(delta / epsilon))
    x = math.pow(2, -1 * J)
    print(f"{delta=}, {epsilon=}, step_size={x}")
    for i in range(0, J + 1):
        depth = int(math.pow(2, i + J))
        print(f"{i=}, {depth=}")
