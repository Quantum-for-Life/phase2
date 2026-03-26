#!/usr/bin/env python3
"""Pauli rotation examples using phase2.

Demonstrates single-qubit, two-qubit, and multi-qubit Pauli
rotations via phase2.run().  Each section derives the expected
result analytically, runs the computation through the C library,
and verifies the output to machine precision.

The function phase2.run() computes

    <psi| prod_k exp(i * delta * coeffs[k] * P_k) |psi>

where P_k are Pauli strings and |psi> is a computational basis
state.  The product is ordered left-to-right: k=0 is applied
first (innermost), k=last is applied last (outermost).

Prerequisites:
    make shared
    pip install .

Run:
    python examples/pauli_rotation.py

With MPI (optional):
    mpirun -n 2 python examples/pauli_rotation.py
"""

import math
import phase2

TOL = 1e-10
passed = 0

# ================================================================
# Section 1: Single-qubit Z rotation
# ================================================================
#
# The Pauli Z matrix is diagonal in the computational basis:
#
#     Z = diag(+1, -1)
#
# with eigenvalues +1 for |0> and -1 for |1>.  Therefore:
#
#     exp(i*a*Z) = diag(exp(ia), exp(-ia))
#                = cos(a)*I + i*sin(a)*Z
#
# Acting on |0>:
#
#     exp(i*a*Z)|0> = exp(ia)|0>
#
# The overlap with |0> is:
#
#     <0|exp(i*a*Z)|0> = exp(ia) = cos(a) + i*sin(a)
#
# We set delta=0.3 and coeffs=[1], so the rotation angle is
# a = delta * coeff = 0.3.

a = 0.3
result = phase2.run(["Z0"], [1.0], a, "0")
expected = complex(math.cos(a), math.sin(a))

print("Section 1: Single-qubit Z rotation")
print(f"  result   = {result}")
print(f"  expected = {expected}")
assert abs(result - expected) < TOL
passed += 1
print("  PASS")
print()

# ================================================================
# Section 2: Single-qubit X rotation
# ================================================================
#
# The Pauli X matrix is the bit-flip operator:
#
#     X = [[0, 1],
#          [1, 0]]
#
# Its matrix exponential is:
#
#     exp(i*a*X) = cos(a)*I + i*sin(a)*X
#
# Acting on |0>:
#
#     exp(i*a*X)|0> = cos(a)|0> + i*sin(a)|1>
#
# The overlap with |0> is:
#
#     <0|exp(i*a*X)|0> = cos(a)
#
# This is purely real because the X rotation mixes |0> and |1>
# symmetrically, and the imaginary part is carried entirely by
# the |1> component.

a = 0.5
result = phase2.run(["X0"], [1.0], a, "0")
expected = complex(math.cos(a), 0.0)

print("Section 2: Single-qubit X rotation")
print(f"  result   = {result}")
print(f"  expected = {expected}")
assert abs(result - expected) < TOL
passed += 1
print("  PASS")
print()

# ================================================================
# Section 3: Two-qubit ZZ rotation
# ================================================================
#
# The tensor product Z tensor Z acts on two qubits.  In the
# computational basis, ZZ is diagonal with eigenvalues:
#
#     |00> -> (+1)(+1) = +1
#     |01> -> (+1)(-1) = -1
#     |10> -> (-1)(+1) = -1
#     |11> -> (-1)(-1) = +1
#
# On the initial state |00>:
#
#     exp(i*a*ZZ)|00> = exp(i*a*1)|00> = exp(ia)|00>
#
# The overlap is identical to the single-qubit Z case:
#
#     <00|exp(i*a*ZZ)|00> = exp(ia)

a = 0.3
result = phase2.run(["Z0 Z1"], [1.0], a, "00")
expected = complex(math.cos(a), math.sin(a))

print("Section 3: Two-qubit ZZ rotation")
print(f"  result   = {result}")
print(f"  expected = {expected}")
assert abs(result - expected) < TOL
passed += 1
print("  PASS")
print()

# ================================================================
# Section 4: Non-commuting X then Z (Trotter step)
# ================================================================
#
# When Pauli operators act on the same qubit, they generally
# do not commute: XZ = -ZX.  The product of exponentials is
# therefore not the exponential of the sum.
#
# phase2.run() applies the operators left to right:
#
#     <0| exp(i*b*Z) exp(i*a*X) |0>
#
# where the first term in the paulis list (X0) is innermost
# and the second (Z0) is outermost.
#
# Step 1: Apply exp(i*a*X) to |0>:
#
#     exp(i*a*X)|0> = cos(a)|0> + i*sin(a)|1>
#
# Step 2: Apply exp(i*b*Z) to the result.  Since Z is
# diagonal, it multiplies each basis state by its eigenvalue:
#
#     exp(i*b*Z)(cos(a)|0> + i*sin(a)|1>)
#       = cos(a)*exp(ib)|0> + i*sin(a)*exp(-ib)|1>
#
# Step 3: Overlap with <0|:
#
#     <0|...> = cos(a) * exp(ib)
#
# With paulis=["X0", "Z0"], coeffs=[1.0, 1.0], delta=0.3:
#   angle_x = delta * 1.0 = 0.3
#   angle_z = delta * 1.0 = 0.3

delta = 0.3
a = delta * 1.0  # X rotation angle
b = delta * 1.0  # Z rotation angle

result = phase2.run(["X0", "Z0"], [1.0, 1.0], delta, "0")
expected = math.cos(a) * complex(math.cos(b), math.sin(b))

print("Section 4: Non-commuting X then Z (Trotter step)")
print(f"  result   = {result}")
print(f"  expected = {expected}")
assert abs(result - expected) < TOL
passed += 1
print("  PASS")
print()

# ================================================================
# Section 5: 10-qubit product of independent X rotations
# ================================================================
#
# When Pauli operators act on distinct qubits, they commute.
# The product of exponentials factorises into a tensor product:
#
#     prod_k exp(i*d*h_k*X_k)
#       = tensor_k exp(i*d*h_k*X_k)
#
# Each factor contributes independently to the overlap:
#
#     <0...0| prod_k exp(i*d*h_k*X_k) |0...0>
#       = prod_k <0|exp(i*d*h_k*X_k)|0>
#       = prod_k cos(d*h_k)
#
# This product is purely real because each X rotation on |0>
# has a purely real overlap (see Section 2).
#
# We use 10 qubits with coefficients h_k = 0.1*(k+1) for
# k = 0, 1, ..., 9, giving h = [0.1, 0.2, ..., 1.0].

nqb = 10
delta = 0.3
h = [0.1 * (k + 1) for k in range(nqb)]

paulis = [f"X{k}" for k in range(nqb)]
psi = "0" * nqb

result = phase2.run(paulis, h, delta, psi)

expected_real = 1.0
for k in range(nqb):
    expected_real *= math.cos(delta * h[k])
expected = complex(expected_real, 0.0)

print("Section 5: 10-qubit product of independent X rotations")
print(f"  result   = {result}")
print(f"  expected = {expected}")
assert abs(result - expected) < TOL
passed += 1
print("  PASS")
print()

# ================================================================
# Section 6: Summary
# ================================================================

print(f"All {passed}/{passed} tests passed.")
