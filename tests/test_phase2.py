"""Extensive pytest suite for the phase2 Python module.

Tests the public API: phase2.run(paulis, coeffs, delta, psi).

The function computes
    <psi| prod_k exp(i * delta * coeffs[k] * P_k) |psi>
where P_k are Pauli strings in "X0 Z1" notation and |psi>
is a computational-basis state given as a bitstring.
"""

import math

import pytest

import phase2

TOL = 1e-12


# -- Return type ------------------------------------------------------

class TestReturnType:
    """Verify that phase2.run returns a Python complex."""

    def test_returns_complex(self):
        result = phase2.run(["Z0"], [1.0], 0.1, "0")
        assert isinstance(result, complex)

    def test_real_part_is_float(self):
        result = phase2.run(["Z0"], [1.0], 0.1, "0")
        assert isinstance(result.real, float)

    def test_imag_part_is_float(self):
        result = phase2.run(["Z0"], [1.0], 0.1, "0")
        assert isinstance(result.imag, float)


# -- Determinism -------------------------------------------------------

class TestDeterminism:
    """Repeated calls with identical arguments must return
    identical results."""

    def test_repeated_calls_equal(self):
        args = (["X0", "Z1"], [0.3, -0.7], 0.25, "01")
        a = phase2.run(*args)
        b = phase2.run(*args)
        assert a == b

    def test_10q_repeated(self):
        paulis = [f"X{k}" for k in range(10)]
        coeffs = [0.1 * (k + 1) for k in range(10)]
        psi = "0" * 10
        a = phase2.run(paulis, coeffs, 0.05, psi)
        b = phase2.run(paulis, coeffs, 0.05, psi)
        assert a == b


# -- Identity (empty Hamiltonian) -------------------------------------

class TestIdentity:
    """An empty product of exponentials is the identity
    operator."""

    def test_empty_paulis(self):
        result = phase2.run([], [], 1.0, "00")
        assert abs(result - 1.0) < TOL

    def test_empty_paulis_single_qubit(self):
        result = phase2.run([], [], 1.0, "0")
        assert abs(result - 1.0) < TOL


# -- Single Pauli exponential -----------------------------------------

class TestSinglePauli:
    """Analytically solvable single-term cases.

    exp(i*a*Z)|0> = (cos a + i sin a)|0>
    exp(i*a*Z)|1> = (cos a - i sin a)|1>
    exp(i*a*X)|0> = cos(a)|0> + i sin(a)|1>
        so <0|exp(i*a*X)|0> = cos(a)
    """

    def test_z_on_0(self):
        a = 0.3
        result = phase2.run(["Z0"], [1.0], a, "0")
        assert abs(result.real - math.cos(a)) < TOL
        assert abs(result.imag - math.sin(a)) < TOL

    def test_z_on_1(self):
        a = 0.7
        result = phase2.run(["Z0"], [1.0], a, "1")
        assert abs(result.real - math.cos(a)) < TOL
        assert abs(result.imag - (-math.sin(a))) < TOL

    def test_x_real(self):
        a = 0.5
        result = phase2.run(["X0"], [1.0], a, "0")
        assert abs(result.real - math.cos(a)) < TOL
        assert abs(result.imag) < TOL

    def test_negative_delta(self):
        a = 0.6
        result = phase2.run(["Z0"], [1.0], -a, "0")
        assert abs(result.real - math.cos(a)) < TOL
        assert abs(result.imag - (-math.sin(a))) < TOL

    def test_delta_scales_coeffs(self):
        r1 = phase2.run(["Z0"], [0.5], 2.0, "0")
        r2 = phase2.run(["Z0"], [1.0], 1.0, "0")
        assert abs(r1 - r2) < TOL


# -- Sequential exponentials -------------------------------------------

class TestSequential:
    """Products of two identical Pauli exponentials."""

    def test_two_x(self):
        """exp(i*a*X) exp(i*a*X) = exp(2iaX), so
        <0|exp(2iaX)|0> = cos(2a)."""
        a = 0.4
        result = phase2.run(
            ["X0", "X0"], [1.0, 1.0], a, "0"
        )
        assert abs(result.real - math.cos(2 * a)) < TOL

    def test_two_z(self):
        """exp(i*a*Z) exp(i*a*Z) = exp(2iaZ), so
        <0|exp(2iaZ)|0> = exp(2ia)."""
        a = 0.35
        result = phase2.run(
            ["Z0", "Z0"], [1.0, 1.0], a, "0"
        )
        expected = complex(
            math.cos(2 * a), math.sin(2 * a)
        )
        assert abs(result - expected) < TOL


# -- Multi-qubit -------------------------------------------------------

class TestMultiQubit:
    """Multi-qubit Pauli strings."""

    def test_zz_on_00(self):
        """Z0 Z1 has eigenvalue +1 on |00>, so
        exp(i*a*Z0Z1)|00> = exp(ia)|00>."""
        a = 0.8
        result = phase2.run(["Z0 Z1"], [1.0], a, "00")
        expected = complex(math.cos(a), math.sin(a))
        assert abs(result - expected) < TOL

    def test_3q_norm(self):
        """The overlap magnitude must not exceed 1."""
        result = phase2.run(
            ["X0 Y1", "Z1 Z2"], [0.3, -0.5], 0.2, "010"
        )
        assert abs(result) <= 1.0 + TOL


# -- Large system (10 qubits) -----------------------------------------

class TestLargeSystem:
    """Ten-qubit tests with analytically known answers."""

    def test_10q_independent_x(self):
        """Ten independent X_k on |0...0>.

        Each factor contributes cos(delta * h_k) to the
        diagonal overlap, so the product is
            prod_k cos(delta * h_k).
        """
        nq = 10
        delta = 0.05
        h = [0.1 * (k + 1) for k in range(nq)]
        paulis = [f"X{k}" for k in range(nq)]
        psi = "0" * nq
        result = phase2.run(paulis, h, delta, psi)
        expected = math.prod(
            math.cos(delta * h[k]) for k in range(nq)
        )
        assert abs(result.real - expected) < TOL
        assert abs(result.imag) < TOL

    def test_10q_x_then_z_per_qubit(self):
        """X_k followed by Z_k on each qubit.

        The product is ordered as
            exp(i*a_0*X_0) exp(i*b_0*Z_0)
            exp(i*a_1*X_1) exp(i*b_1*Z_1) ...

        For qubit k in state |0>:
            <0| exp(i*a_k*X_k) exp(i*b_k*Z_k) |0>
          = cos(a_k) * exp(i*b_k)

        The total overlap is
            prod_k cos(a_k) * exp(i * sum_k b_k).
        """
        nq = 10
        a = [0.01 * (k + 1) for k in range(nq)]
        b = [0.02 * (k + 1) for k in range(nq)]
        paulis = []
        coeffs = []
        for k in range(nq):
            paulis.append(f"X{k}")
            coeffs.append(a[k])
            paulis.append(f"Z{k}")
            coeffs.append(b[k])
        psi = "0" * nq
        result = phase2.run(paulis, coeffs, 1.0, psi)
        mag = math.prod(math.cos(a[k]) for k in range(nq))
        phase = sum(b[k] for k in range(nq))
        expected = mag * complex(
            math.cos(phase), math.sin(phase)
        )
        assert abs(result - expected) < TOL

    def test_10q_z_on_alternating(self):
        """Z_k on |0101010101>.

        Z has eigenvalue +1 on |0> and -1 on |1>.
        For the alternating state the eigenvalue of Z_k is
            sign_k = +1 if bit k is 0, -1 if bit k is 1.
        The overlap is
            exp(i * sum_k sign_k * delta * h_k).
        """
        nq = 10
        delta = 0.1
        h = [0.1 * (k + 1) for k in range(nq)]
        psi = "0101010101"
        paulis = [f"Z{k}" for k in range(nq)]
        signs = [1 if psi[k] == "0" else -1
                 for k in range(nq)]
        result = phase2.run(paulis, h, delta, psi)
        phase = sum(
            signs[k] * delta * h[k] for k in range(nq)
        )
        expected = complex(
            math.cos(phase), math.sin(phase)
        )
        assert abs(result - expected) < TOL


# -- Error handling ----------------------------------------------------

class TestErrors:
    """Invalid inputs must raise ValueError."""

    def test_length_mismatch(self):
        with pytest.raises(ValueError):
            phase2.run(["Z0"], [1.0, 2.0], 0.1, "0")

    def test_invalid_pauli_operator(self):
        with pytest.raises(ValueError):
            phase2.run(["W0"], [1.0], 0.1, "0")

    def test_qubit_out_of_range(self):
        with pytest.raises(ValueError):
            phase2.run(["X5"], [1.0], 0.1, "0")

    def test_invalid_bitstring(self):
        with pytest.raises(ValueError):
            phase2.run(["Z0"], [1.0], 0.1, "2")


# -- Edge cases --------------------------------------------------------

class TestEdgeCases:
    """Boundary values of delta."""

    def test_zero_delta(self):
        result = phase2.run(["X0"], [1.0], 0.0, "0")
        assert abs(result - 1.0) < TOL

    def test_delta_pi(self):
        """exp(i*pi*Z)|0> = (cos pi + i sin pi)|0>
        = -1."""
        result = phase2.run(
            ["Z0"], [1.0], math.pi, "0"
        )
        assert abs(result - (-1.0)) < TOL

    def test_delta_half_pi(self):
        """exp(i*pi/2*Z)|0> = (cos pi/2 + i sin pi/2)|0>
        = i."""
        result = phase2.run(
            ["Z0"], [1.0], math.pi / 2, "0"
        )
        assert abs(result - 1j) < TOL


# -- Cross-validation --------------------------------------------------

