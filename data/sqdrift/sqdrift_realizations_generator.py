"""Minimal SqDRIFT mock: molecule -> qDRIFT realizations -> list of qiskit circuits.

Pure qiskit-fermions + qiskit + pyscf. Steps:
  1. PySCF CASCI active space -> fermionic Hamiltonian.
  2. Circuit = Hartree-Fock state + one exact evolution gate exp(-i H t).
  3. qDRIFT replaces that gate by N_EXCITATIONS randomly sampled ones,
     once per realization (each with its own seed).
  4. Jordan-Wigner maps every realization to a QuantumCircuit.

Defaults: N2 cc-pVDZ CAS(6e,6o) -> 12 qubits. Edit the constants below.

Run single-threaded on a login node -- the CPU-time cap adds up over threads,
so PySCF's default thread fan-out trips it within seconds:

    OMP_NUM_THREADS=1 python mock_sqdrift_qiskit.py
"""

from __future__ import annotations

import numpy as np
from pyscf import ao2mo, gto, lib, mcscf, scf
from qiskit import QuantumCircuit
from qiskit_fermions.circuit import FermionicCircuit
from qiskit_fermions.circuit.library import Evolution, InitializeModes
from qiskit_fermions.operators import FermionOperator
from qiskit_fermions.operators.terms.grouping import group_terms_by_electronic_structure
from qiskit_fermions.transpiler.converters import FermionicCircuitToDAG, FermionicDAGToCircuit
from qiskit_fermions.transpiler.passes import QDriftTrotterization
from qiskit_fermions.transpiler.presets import generate_preset_jw_pass_manager

# General settings
BASIS = "cc-pvdz"
NCAS = 6
NELEC = (3, 3)
MO_WINDOW = [4, 5, 6, 7, 8, 16]  # hand-picked N2 CAS(6e,6o) orbitals, 0-based

EVOLUTION_TIME = 1.0 # (t)
N_EXCITATIONS = 10  # (N_E) qDRIFT samples per realization
N_REALIZATIONS = 3  # (N_R) instances generated in this run

OPTIMIZATION_LEVEL = 1
BASIS_GATES = ["cx", "rz", "sx", "x"]


if __name__ == "__main__":
    num_modes = 2 * NCAS

    # active space hamiltonian
    print(f"[1] N2/{BASIS} CAS({sum(NELEC)}e,{NCAS}o) -> {num_modes} qubits", flush=True)
    mol = gto.M(atom="N 0 0 0; N 0 0 1.098", basis=BASIS, symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    mycas = mcscf.CASCI(mf, NCAS, NELEC)
    e_casci = mycas.kernel(mycas.sort_mo(MO_WINDOW, base=0))[0]
    print(f"    RHF   energy: {mf.e_tot:.6f} Ha")
    print(f"    CASCI energy: {e_casci:.6f} Ha", flush=True)

    h1e, ecore = mycas.get_h1eff()
    h2e = mycas.get_h2eff()

    # qiskit-fermions builds the operator straight from the packed integrals.
    # ecore is left out on purpose: it is only a global phase, but as an identity
    # term qDRIFT could sample it and waste a gate on nothing.
    hamiltonian = FermionOperator.from_1body_tril_spin_sym(lib.pack_tril(h1e), norb=NCAS)
    hamiltonian += FermionOperator.from_2body_tril_spin_sym(
        ao2mo.restore(8, np.ascontiguousarray(h2e), NCAS), norb=NCAS
    )

    # Normal-order and drop the identity term. Sorting keeps the term order (and
    # so the sampling) reproducible: simplify() does not guarantee an order.
    hamiltonian = hamiltonian.normal_ordered().simplify(atol=1e-16)
    hamiltonian = FermionOperator.from_terms(
        sorted(
            ((actions, coeff) for actions, coeff in hamiltonian.iter_terms() if actions),
            key=lambda tc: tc[0],
        )
    )
    # Pairs each term with its Hermitian conjugate (in place, returns None). These
    # groups are what qDRIFT samples, so every sampled gate conserves electrons.
    assert (
        group_terms_by_electronic_structure(
            hamiltonian, num_modes, two_body_physicist_order=False
        )
        is None
    )
    print(f"    terms (identity dropped):  {len(hamiltonian)}")
    print(f"    distinct Hermitian groups: {len(set(hamiltonian.groups))}", flush=True)

    # FermionicCircuit
    print(f"[2] FermionicCircuit: HF{NELEC} + Evolution(H, t={EVOLUTION_TIME})", flush=True)
    circ = FermionicCircuit(num_modes)
    # HF reference: lowest n_alpha modes of the alpha sector, lowest n_beta of the
    # beta sector, in block-spin order (alpha 0..NCAS-1, beta NCAS..).
    circ.append(InitializeModes.from_hartree_fock(NCAS, NELEC), circ.modes)
    circ.append(Evolution(num_modes, hamiltonian, EVOLUTION_TIME), circ.modes)
    # No randomness in here yet, so this is built once and re-sampled below.
    dag = FermionicCircuitToDAG().run(circ)

    # One realization per seed
    print(
        f"[3] qDRIFT: {N_REALIZATIONS} realizations x {N_EXCITATIONS} excitations",
        flush=True,
    )
    # Plain Jordan-Wigner: no device coupling map, no mode reordering.
    pm = generate_preset_jw_pass_manager(
        optimization_level=OPTIMIZATION_LEVEL, basis_gates=BASIS_GATES
    )

    circuits: list[QuantumCircuit] = []
    header = f"{'idx':>4} | {'evolutions':>10} | {'qubits':>6} | {'depth':>6} | {'2q gates':>8}"
    print(header)
    print("-" * len(header))
    for idx in range(N_REALIZATIONS):
        # One seeded generator per realization, so any single realization can be
        # reproduced on its own (e.g. in a SLURM array task) without replaying
        # the ones before it. filter_diagonal_terms drops number-operator terms:
        # they only add phases, which the measured bitstrings never see.
        sampled_dag = QDriftTrotterization(
            N_EXCITATIONS, filter_diagonal_terms=True, rng=idx
        ).run(dag)
        n_evo = sum(1 for node in sampled_dag.op_nodes() if isinstance(node.op, Evolution))

        # The pass manager expects a FermionicCircuit, not a DAG.
        qc = pm.run(FermionicDAGToCircuit().run(sampled_dag))
        circuits.append(qc)

        print(
            f"{idx:>4} | {n_evo:>10} | {qc.num_qubits:>6} | {qc.depth():>6} | "
            f"{qc.num_nonlocal_gates():>8}",
            flush=True,
        )

    print(f"\n[4] {len(circuits)} QuantumCircuits generated.")
    print(f"    realization 0 ops: {dict(circuits[0].count_ops())}")
