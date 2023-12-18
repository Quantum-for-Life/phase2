import warnings

warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
import numpy as np
from pyscf import gto, mcscf
from pyscf.tools import fcidump
from qiskit_nature.second_q.formats import fcidump_to_problem
from qiskit_nature.second_q.formats.fcidump import FCIDump
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library.initial_states.hartree_fock import (
    hartree_fock_bitstring,
)
import h5py

PAULI_DICT = {
    'I': 0,
    'X': 1,
    'Y': 2,
    'Z': 3,
}

def main():
    xyz = """O   0.000000    0.000000    0.117790
             H   0.000000    0.755453    -0.471161
             H   0.000000    -0.755453   -0.471161"""

    # scf calculation
    mol = gto.M(atom=xyz)
    mol.basis = "sto3g"
    mol.spin = 0
    mol.build()
    hf = mol.RHF().run()
    # print(f"\nSCF energy: {hf.e_tot}")
    # print(hf.analyze())
    e_nuc = hf.energy_nuc()  # nuclear repulsion energy

    filename = "H2O_CAS56"
    fci_file = filename + ".FCI_DUMP"
    h5file = filename + ".h5"
    cas = (5, 6)
    casci = mcscf.CASCI(hf, cas[0], cas[1])
    res_casci = casci.kernel()
    e_tot = res_casci[0]  # el. energy in CAS + el. energy in core + nuclear repulsion energy
    e_cas = res_casci[1]  # el. energy in CAS
    e_const = e_tot - e_cas  # el. energy in core + nuclear repulsion energy

    fcidump.from_mcscf(casci, fci_file)

    fcidump_qiskit = FCIDump.from_file(fci_file)
    el_str_problem = fcidump_to_problem(fcidump_qiskit)

    hamiltonian = el_str_problem.hamiltonian
    fermionic_op = hamiltonian.second_q_op()

    mapper = JordanWignerMapper()
    qubit_op = mapper.map(fermionic_op)

    # eigenvalues of Hamiltonian after JW
    H_matrix = np.real(qubit_op.to_matrix())
    eigvals = np.sort(np.linalg.eigvals(H_matrix))
    # print("diagonalization of qubit Hamiltonian:")
    # pprint.pprint(set(eigvals))
    TOL = 1e-6
    assert abs(min(eigvals) - res_casci[1]) < TOL
    print("Qubit Hamiltonian preserved lowest eigenvalue. Success!")

    # HF state
    half_elec = int(cas[1] / 2)
    assert 2 * half_elec == cas[1]  # assert even number of active electrons
    hf_bitstring = hartree_fock_bitstring(cas[0], (half_elec, half_elec))
    hf_bitstring = [1 if bit else 0 for bit in hf_bitstring]

    # preparing data to write to h5 file
    num_qubits = qubit_op.num_qubits
    qubit_op_list = qubit_op.to_list()
    offset = 0
    identity_index = None
    for i,(label, coeff) in enumerate(qubit_op.to_list()):
        if label == num_qubits * "I":
            offset = coeff.real
            identity_index = i
            break

    if identity_index is None:
        pass
    else:
        del qubit_op_list[identity_index]
    num_terms = len(qubit_op_list)

    pauli_labels = np.array([term[0] for term in qubit_op_list])
    pauli_coeffs = np.array([term[1] for term in qubit_op_list])
    assert np.allclose(np.imag(pauli_coeffs), 0)
    pauli_coeffs = np.real(pauli_coeffs)
    pauli_labels_enc = []
    for label in pauli_labels:
        label_enc = []
        for c in label[::-1]:
            label_enc.append(PAULI_DICT[c])
        pauli_labels_enc.append(label_enc)
    normalization = 1/np.sum(np.abs(pauli_coeffs))

    # write to h5 file
    with h5py.File(h5file, "w") as f:
        h5_ferm_bm = f.create_group("fermion/benchmark")
        h5_ferm_bm.attrs["method"] = "CASCI(5,6)"
        h5_ferm_bm.attrs["energy"] = e_tot

        h5_ph = f.create_group("pauli_hamil")
        h5_ph.create_dataset("coeffs", shape=(num_terms,), dtype="d")[...] = pauli_coeffs
        h5_ph.create_dataset("paulis", shape=(num_terms, num_qubits), dtype="u1")[...] = pauli_labels_enc
        h5_ph.attrs["normalization"] = normalization
        h5_ph.attrs["offset"] = offset + e_const

        h5_sp = f.create_group("state_prep/multidet")
        h5_sp.create_dataset("coeffs", shape=(1,2), dtype="d")[...] = [[1, 0]]
        h5_sp.create_dataset("dets", shape=(1, num_qubits), dtype="u1")[...] = [hf_bitstring]

if __name__ == "__main__":
    main()
