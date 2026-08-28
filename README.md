# Molecular Quantum Simulations for Drug Metabolism Optimization

Data accompanying the **Open Quantum Institute (OQI) Phase 3 report** of the
Algorithmiq team.

Cytochrome P450 (CYP) enzymes govern the metabolism of a large fraction of
clinically used drugs. The rate of a metabolic step depends exponentially on its
activation barrier, so predicting metabolism quantitatively requires
electronic-structure accuracy that is hard to reach for the open-shell iron-oxo
chemistry of the reactive species Compound I (Cpd I).

This project benchmarks a hybrid quantum-classical workflow on a controlled
model of that chemistry: the **aromatic hydroxylation of benzene by Cpd I**. The
target quantity is the activation barrier

```
ΔE‡ = E(TS) − E(RC)
```

between the reactant complex (RC) and the transition state (TS). Benzene stands
in for pharmaceutically relevant aromatic substrates; the intended follow-on
application is S-warfarin metabolised by CYP2C9.

> ### Disclaimer
>
> This repository provides the **intermediate results required to repeat our
> calculations**: geometries, converged mean-field solutions, active-space
> Hamiltonians and quantum circuits. Results and workflow descriptions are added
> **progressively** as the project advances, so the contents here are at any
> given time a partial view of the study.
>
> **The Phase 3 report remains the authoritative and up-to-date overview of all
> relevant results.**

## Contents

- [1. Data](#1-data)
  - [1.1 Geometries](#11-geometries)
  - [1.2 UHF checkpoints](#12-uhf-checkpoints)
  - [1.3 Active-space Hamiltonians](#13-active-space-hamiltonians)
  - [1.4 ADAPT-VMPE circuits](#14-adapt-vmpe-circuits)
- [2. Methods](#2-methods)
  - [2.1 Classical characterisation of the reaction pathway](#21-classical-characterisation-of-the-reaction-pathway)
  - [2.2 UHF reference state](#22-uhf-reference-state)
  - [2.3 FNO active spaces and Hamiltonian construction](#23-fno-active-spaces-and-hamiltonian-construction)
  - [2.4 CCSD benchmarking of the active spaces](#24-ccsd-benchmarking-of-the-active-spaces)
  - [2.5 DMRG reference calculations](#25-dmrg-reference-calculations)
  - [2.6 ADAPT-VMPE circuit construction](#26-adapt-vmpe-circuit-construction)
  - [2.7 SqDRIFT sampling and QSCI](#27-sqdrift-sampling-and-qsci)
- [3. Results](#3-results)
  - [3.1 DMRG](#31-dmrg)

The chain of dependencies is linear, with the active-space Hamiltonian as its
central node:

```
  geometries ──▶ UHF reference ──▶ FNO active space ──▶ Hamiltonian ──┬──▶ DMRG reference
                                                                      ├──▶ ADAPT-VMPE circuits
                                                                      └──▶ SqDRIFT → QSCI
```

# 1. Data

## 1.1 Geometries

`data/geometries/`

The structures are taken as published from the Supporting Information of

> Lonsdale, R.; Harvey, J. N.; Mulholland, A. J.
> *Effects of Dispersion in Density Functional Based Quantum
> Mechanical/Molecular Mechanical Calculations on Cytochrome P450 Catalyzed
> Reactions.* J. Chem. Theory Comput. **8**, 4637 (2012).

That study treated the addition of benzene to Cpd I in CYP2C9 with
dispersion-corrected QM/MM. No re-optimisation was performed here.

The model comprises the Fe=O unit, the porphyrin macrocycle, an axial SCH₃ group
representing the cysteine thiolate ligand, and the benzene substrate: 55 atoms,
neutral overall. Eight structures cover all four reaction pathways at both ends
of the reaction coordinate:

| | face-on doublet | face-on quartet | side-on doublet | side-on quartet |
|---|---|---|---|---|
| **RC** | `RC_face_on_doublet.xyz` | `RC_face_on_quartet.xyz` | `RC_side_on_doublet.xyz` | `RC_side_on_quartet.xyz` |
| **TS** | `TS_face_on_doublet.xyz` | `TS_face_on_quartet.xyz` | `TS_side_on_doublet.xyz` | `TS_side_on_quartet.xyz` |

The two orientations differ in how benzene approaches the heme: *side-on* has
the aromatic ring approximately perpendicular to the porphyrin plane, *face-on*
approximately parallel. The doublet and quartet are the two low-lying spin
manifolds of Cpd I, which carries three open-shell electrons.

Files are standard XYZ in Ångström, with charge and multiplicity in the comment
line.

### Why the side-on doublet

All correlated and quantum-computing work uses the **side-on doublet** pathway.
That choice follows from the benchmark in
[Section 2.1](#21-classical-characterisation-of-the-reaction-pathway), where the
barrier of all four pathways was computed with HF, MP2, B3LYP and B3LYP-D4
across several basis sets.

HF and MP2 barriers turned out to depend strongly on both spin state and basis
set, as expected for a system with several coupled open-shell orbitals. The DFT
barriers were far more stable, and among them the side-on doublet gave the
lowest barrier, near 14 kcal/mol for B3LYP-D4. Since the reaction rate depends
exponentially on the barrier, the lowest accessible pathway is expected to
dominate product formation.

## 1.2 UHF checkpoints

`data/uhf_checkpoints/`

| File | Content |
|---|---|
| `RC_side_on_doublet_uhf.chk` | converged UHF/cc-pVDZ orbitals of the RC, PySCF checkpoint |
| `input_rc.py` | reaches the RC UHF state from that checkpoint |
| `input_ts.py` | reaches the TS UHF state from that checkpoint |

The doublet ground state of Cpd I is a broken-symmetry solution, a triplet
iron-oxo unit coupled antiferromagnetically to a ligand radical. An SCF started
from a default guess is not guaranteed to converge to it. The same checkpoint
seeds both the RC and the TS, so that both converge to the same electronic
solution and their energy difference is a meaningful barrier.

```bash
cd data/uhf_checkpoints
python input_rc.py     # or input_ts.py
```

Each script reads its geometry from `data/geometries/`, builds the initial
density matrix from the checkpoint, and converges a second-order (Newton) UHF.
The scripts write their own checkpoints and leave the published file untouched.

To inspect the stored orbitals without running an SCF:

```python
from pyscf import lib, scf

chk = "RC_side_on_doublet_uhf.chk"
mol = lib.chkfile.load_mol(chk)          # the molecule travels with the checkpoint
data = scf.chkfile.load(chk, "scf")      # e_tot, mo_coeff, mo_energy, mo_occ

mo_a, mo_b = data["mo_coeff"]            # unrestricted: one set per spin
```

⟨S²⟩ of the converged solution lies well above the spin-pure doublet value of
0.75, around 1.75 in the idealised three-open-shell limit. This reflects the
antiferromagnetically coupled three-spin structure of Cpd I rather than ordinary
spin contamination.

## 1.3 Active-space Hamiltonians

`data/as_hamiltonians/`

Second-quantised Hamiltonians for the side-on doublet RC and TS, in a hierarchy
of Frozen Natural Orbital (FNO) active spaces:

```
H = E_core + Σ_pq t_pq a†_p a_q + ½ Σ_pqrs v_pqrs a†_p a†_q a_s a_r
```

`E_core` contains the nuclear repulsion and the frozen-orbital contribution, so
total energies are directly comparable to the reference values. RC and TS have
different `E_core`, so only **total** energies may be subtracted to form the
barrier.

| Directory | Active space | Spatial orbitals | Electrons (α, β) | Qubits | `ham.npz` |
|---|---|---|---|---|---|
| `RC_11_10`, `TS_11_10` | CAS(11,10) | 10 | (6, 5) | 20 | 0.2 MB |
| `RC_21_20`, `TS_21_20` | CAS(21,20) | 20 | (11, 10) | 40 | 3.8 MB |
| `RC_31_30`, `TS_31_30` | CAS(31,30) | 30 | (16, 15) | 60 | 20 MB |
| `RC_41_40`, `TS_41_40` | CAS(41,40) | 40 | (21, 20) | 80 | 62 MB |
| *(not on GitHub)* | CAS(51,50) | 50 | (26, 25) | 100 | 150 MB |

The CAS(51,50) integrals exceed the 100 MiB per-file limit that GitHub imposes
and are therefore not hosted here. They are available on request; please open an
issue or contact the authors.

CAS(11,10) is the benchmark case: its symmetry sector holds 52 920 determinants,
small enough for an exact unrestricted CASCI (UCASCI) reference.

Each directory holds two files:

| File | Key | Shape | Content |
|---|---|---|---|
| `ham.npz` | `h0` | scalar | core energy `E_core` / Hartree |
| | `h1` | `(2, n, n)` | one-electron integrals, ordered `(α, β)` |
| | `h2` | `(3, n, n, n, n)` | two-electron integrals, ordered `(αα, αβ, ββ)` |
| `overlap.npz` | `mo_alpha` | `(n_ao, n)` | active α MO coefficients, `n_ao = 614` |
| | `mo_beta` | `(n_ao, n)` | active β MO coefficients |
| | `ao_overlap` | `(n_ao, n_ao)` | AO overlap matrix |

Reading a Hamiltonian:

```python
import numpy as np

ham = np.load("data/as_hamiltonians/RC_11_10/ham.npz")

e_core = ham["h0"].item()      # scalar
h1     = ham["h1"]             # (2, n, n)        -> alpha, beta
h2     = ham["h2"]             # (3, n, n, n, n)  -> aa, ab, bb

n_orb = h1.shape[-1]
print(f"{n_orb} spatial orbitals, E_core = {e_core:.9f} Ha")
```

The integrals follow the convention of `pyblock2`'s `get_uhf_integrals`, so they
can be passed to a block2 driver directly:

```python
from pyblock2.driver.core import DMRGDriver, SymmetryTypes

driver = DMRGDriver(scratch="./tmp", symm_type=SymmetryTypes.SZ, n_threads=8)
driver.initialize_system(n_sites=n_orb, n_elec=11, spin=1)   # CAS(11,10), 2Sz = 1
mpo = driver.get_qc_mpo(h1e=h1, g2e=h2, ecore=e_core, integral_cutoff=1e-12)
```

`overlap.npz` is provided separately because the active orbitals are UHF
orbitals, whose α and β sets are not mutually orthogonal. Quantities depending
on that overlap, ⟨S²⟩ above all, cannot be evaluated from the integrals alone:

```python
from pyscf.fci.spin_op import spin_square_general

ov = np.load("data/as_hamiltonians/RC_11_10/overlap.npz")
ssq, mult = spin_square_general(
    dm1a, dm1b, dm2aa, dm2ab, dm2bb,
    (ov["mo_alpha"], ov["mo_beta"]), ov["ao_overlap"],
)
```

## 1.4 ADAPT-VMPE circuits

`data/adapt_circuits/`

ADAPT-VMPE circuits for the CAS(11,10) RC and TS (20 qubits), as
`rc_11e20q.tar.gz` and `ts_11e20q.tar.gz`.

*Documentation of the file format and the larger active spaces will be added.*

# 2. Methods

How each calculation is carried out, starting from the data published here. The
core solvers (the in-house distributed multi-GPU DMRG implementation,
ADAPT-VMPE, SqDRIFT and Treespilation) are proprietary and are not part of this
release. The published data are nevertheless sufficient to check the scientific
claims: the Hamiltonians are self-contained, and any external solver can be
benchmarked against the reference energies.

## 2.1 Classical characterisation of the reaction pathway

*Input: `data/geometries/`, all eight structures.*

- Run single-point calculations on the RC and TS of each pathway with HF, MP2,
  B3LYP and B3LYP-D4, across a series of basis sets.
- Form the barrier `E_a = E(TS) − E(RC)` for every combination of method, basis
  set, orientation and spin state.
- Compute Löwdin spin populations at the B3LYP/cc-pVDZ level for the side-on
  doublet, grouped into the fragments FeO, SCH₃, porphyrin and benzene, to
  follow how spin is redistributed as the C-O bond forms.
- Select the pathway with the lowest DFT barrier for all further work; see
  [Why the side-on doublet](#why-the-side-on-doublet).

## 2.2 UHF reference state

*Input: `data/geometries/{RC,TS}_side_on_doublet.xyz`,
`data/uhf_checkpoints/RC_side_on_doublet_uhf.chk`.*

- Build the molecule at charge 0 and `spin = 1` in the cc-pVDZ basis.
- Take the initial density matrix from the published checkpoint, using the same
  checkpoint for RC and TS.
- Converge a second-order (Newton) UHF to `conv_tol = 1e-8`.
- Check ⟨S²⟩ against the broken-symmetry reference of ≈ 1.75, not against the
  spin-pure 0.75.

This is what `input_rc.py` and `input_ts.py` do; run them directly.

## 2.3 FNO active spaces and Hamiltonian construction

*Input: the converged UHF solution. Output: `data/as_hamiltonians/`.*

- Select the occupied active orbitals by orbital energy, keeping those closest
  to the Fermi level.
- Compress the virtual space with the FNO procedure: run MP2 on the canonical HF
  orbitals, build and diagonalise the virtual-virtual block of the MP2
  one-particle density matrix, and keep the natural orbitals with the largest
  occupation numbers.
- Repeat for 10, 20, 30, 40 and 50 active orbitals, giving CAS(11,10) through
  CAS(51,50).
- Transform the integrals into the active MO basis and fold the nuclear
  repulsion and frozen-orbital contribution into `E_core`. Store `h0`, `h1`,
  `h2` alongside the active MO coefficients and the AO overlap.

## 2.4 CCSD benchmarking of the active spaces

*Input: the UHF solution and the FNO active spaces.*

- Compute full-space MP2/cc-pVDZ and CCSD/cc-pVDZ barriers as reference lines
  (density fitting, frozen core), and record the CCSD `T₁` diagnostic.
- For every active space, run three calculations:
  - **FNO-MP2**, MP2 inside the truncated space;
  - **FNO-CCSD**, CCSD in the same space;
  - **FNO-CCSD + MP2corr**, FNO-CCSD plus an MP2 estimate of the correlation
    excluded by the truncation,
    `E = E(FNO-CCSD) + [E(MP2, full) − E(FNO-MP2)]`.
- Compare the barrier, not only the absolute energies. The uncorrected barriers
  do not converge monotonically with active-space size, because the correlation
  omitted by truncation does not cancel between RC and TS.

## 2.5 DMRG reference calculations

*Input: `data/as_hamiltonians/{RC,TS}_*/ham.npz` and `overlap.npz`.*

- Load `h0`, `h1`, `h2` and build the quantum-chemistry MPO; the integrals can
  be passed to a block2 `DMRGDriver` in `SymmetryTypes.SZ` directly.
- Target the ground state at fixed electron number and spin projection
  `2Sz = 1`.
- Converge a ladder of increasing bond dimensions for each state and active
  space, from D = 200 for CAS(11,10) up to D = 7000 for CAS(51,50).
- Extrapolate the energy linearly against the discarded weight δ to the FCI
  limit δ → 0.
- Evaluate ⟨S²⟩ with `spin_square_general`, passing the MO coefficients and AO
  overlap from `overlap.npz`. The naive site-spin operator gives the wrong
  answer, since the α and β active orbitals are not mutually orthogonal.
- A Fiedler orbital reordering, computed from the spin-summed exchange
  integrals, improves convergence at a given bond dimension.

## 2.6 ADAPT-VMPE circuit construction

*Input: `data/as_hamiltonians/`.*

- Start from the Hartree-Fock determinant, augmented by active orbital rotations
  that are re-optimised alongside every operator added.
- Grow the ansatz one operator at a time in the Schrödinger picture: score the
  pool, append the element that most improves the variational energy, and
  re-optimise all parameters afterwards.
- Evaluate energies and gradients with Majorana Propagation at a monomial-length
  cutoff of 6, which keeps the cost polynomial in the number of orbitals.
- Trim the pool: score the complete pool once every 20 iterations, keep the top
  5 %, and evaluate only that subset in between.
- Report errors against the extrapolated DMRG references. Track the barrier
  error, not only the absolute errors, since much of the systematic bias cancels
  in the difference.

## 2.7 SqDRIFT sampling and QSCI

*Input: `data/as_hamiltonians/{RC,TS}_11_10/`, the 20-qubit case with an exact
UCASCI reference, so every reported error is purely a subspace-truncation
error.*

- Partition the normal-ordered Hamiltonian in the fermionic excitation basis and
  remove the diagonal number-operator terms from the sampling pool. Since every
  drawn term is a genuine fermionic excitation, each circuit conserves particle
  number and never leaves the target sector.
- Draw qDRIFT products of `N_E ∈ {5, 10, 25, 50, 100}` random excitations at
  evolution time `t = 1`, with `N_R = 50` seeded realisations per `N_E` and
  `N_S = 1024` shots each.
- Map the fermionic ansatz to hardware and score each pipeline by two-qubit gate
  count and circuit depth. Five were compared: plain Jordan-Wigner with Qiskit
  routing, Jordan-Wigner on a path graph, two simulated-annealing mode
  optimisations, and full Treespilation over Bonsai trees.
- Sample with Qiskit Aer's matrix-product-state simulator at bond dimension
  χ = 512, under a noise model built from a calibration snapshot of IBM's
  156-qubit Heron r2 device. Run once without noise to certify that the counts
  are physical and the χ sweep has converged.
- Filter the sampled bitstrings to the sector `(N_α, N_β) = (6, 5)`, deduplicate,
  and diagonalise `H` in the resulting subspace by matrix-free Lanczos over
  sparse Slater-Condon matrix elements.
- For the dressed variant, fold the orbital rotation of an ADAPT-VMPE ansatz
  into the integrals before sampling and sweep how many gates are dressed.
  Dressing is a one-particle basis change and leaves the exact spectrum
  unchanged, but rotates the orbitals so that the sampled subspace becomes more
  informative at fixed `N_E`.

# 3. Results

## 3.1 DFT, HF, and MP2

Classical HF, MP2, and DFT calculations were performed for the side-on
and face-on reaction pathways in both the doublet and quartet spin states.

### Side-on doublet

| Method | RC (Hartree) | TS (Hartree) | ΔE (kcal/mol) | <S²> RC | <S²> TS |
|---|---:|---:|---:|---:|---:|
| HF/def2-SVP TS from RC guess | -2985.90124 | -2985.89477 | 4.06 | 5.81 | 6.17 |
| HF/def2-SVP RIJCOSX def2/J guess from NoRI | -2985.90256 | -2985.89599 | 4.12 | 5.81 | 6.17 |
| RI-MP2/def2-SVP | -2990.49519 | -2990.47924 | 10.01 | | |
| HF/def2-TZVP | -2987.64656 | -2987.63656 | 6.27 | 5.76 | 6.18 |
| HF/def2-TZVP RIJCOSX def2/J from B3LYP guess | -2987.64764 | -2987.63773 | 6.22 | 5.76 | 6.18 |
| RI-MP2/def2-TZVP | -2993.32687 | -2993.30133 | 16.03 | | |
| HF/def2-QZVP | -2987.74907 | -2987.73830 | 6.76 | 5.77 | 6.20 |
| HF/def2-QZVP RIJCOSX def2/J | -2987.75027 | -2987.73948 | 6.77 | 5.77 | 6.20 |
| RI-MP2/def2-QZVP | -2993.95751 | -2993.93098 | 16.65 | | |
| HF/cc-pVDZ | -2987.32273 | -2987.31440 | 5.23 | 5.78 | 6.16 |
| HF/cc-pVDZ RIJCOSX def2/J | -2987.32220 | -2987.31558 | 4.15 | 5.75 | 6.16 |
| RI-MP2/cc-pVDZ | -2992.05989 | -2992.03934 | 12.89 | | |
| HF/cc-pVTZ | -2987.64780 | -2987.63924 | 5.37 | 5.73 | 6.20 |
| HF/cc-pVTZ RIJCOSX def2/J | -2987.64917 | -2987.64047 | 5.46 | 5.73 | 6.20 |
| RI-MP2/cc-pVTZ | -2993.52549 | -2993.49853 | 16.92 | | |
| HF/cc-pVQZ | -2987.73473 | -2987.72582 | 5.59 | 5.72 | 6.20 |
| HF/cc-pVQZ RIJCOSX def2/J | -2987.73587 | -2987.72688 | 5.65 | 5.72 | 6.20 |
| RI-MP2/cc-pVQZ | -2994.04154 | -2994.01315 | 17.82 | | |
| B3LYP/def2-SVP | -2995.46914 | -2995.43998 | 18.30 | 1.77 | 1.73 |
| B3LYP/def2-TZVP | -2997.21987 | -2997.19120 | 17.99 | 1.77 | 1.71 |
| B3LYP/def2-QZVP | -2997.34912 | -2997.31968 | 18.48 | 1.77 | 1.71 |
| B3LYP/def2-SVP RIJCOSX def2/J | -2995.47048 | -2995.44127 | 18.33 | 1.77 | 1.73 |
| B3LYP/def2-TZVP RIJCOSX def2/J | -2997.22097 | -2997.19232 | 17.98 | 1.77 | 1.71 |
| B3LYP/def2-QZVP RIJCOSX def2/J | -2997.34912 | -2997.32104 | 17.62 | 1.77 | 1.71 |
| B3LYP/cc-pVDZ | -2996.84093 | -2996.81317 | 17.42 | 1.77 | 1.71 |
| B3LYP/cc-pVTZ | -2997.22282 | -2997.19460 | 17.71 | 1.77 | 1.71 |
| B3LYP/cc-pVQZ | -2997.32537 | -2997.29739 | 17.56 | 1.77 | 1.71 |
| B3LYP/cc-pVDZ RIJCOSX def2/J | -2996.84211 | -2996.81437 | 17.41 | 1.77 | 1.71 |
| B3LYP/cc-pVTZ RIJCOSX def2/J | -2997.22415 | -2997.19586 | 17.75 | 1.77 | 1.71 |
| B3LYP/cc-pVQZ RIJCOSX def2/J | -2997.32664 | -2997.29875 | 17.51 | 1.77 | 1.71 |
| B3LYP-D4/def2-SVP RIJCOSX def2/J | -2995.63485 | -2995.61179 | 14.47 | 1.77 | 1.73 |
| B3LYP-D4/def2-TZVP RIJCOSX def2/J | -2997.38534 | -2997.36284 | 14.12 | 1.77 | 1.71 |
| B3LYP-D4/def2-QZVP RIJCOSX def2/J | -2997.51349 | -2997.49157 | 13.76 | 1.77 | 1.71 |
| B3LYP-D4/cc-pVDZ RIJCOSX def2/J | -2997.00648 | -2996.98489 | 13.55 | 1.77 | 1.71 |
| B3LYP-D4/cc-pVTZ RIJCOSX def2/J | -2997.38852 | -2997.36639 | 13.89 | 1.77 | 1.71 |
| B3LYP-D4/cc-pVQZ RIJCOSX def2/J | -2997.49101 | -2997.46927 | 13.64 | 1.77 | 1.71 |
| B3LYP/G/cc-pVDZ | -2997.81448 | -2997.78681 | 17.36 | 1.77 | 1.70 |
| B3LYP/G/cc-pVTZ | -2998.19625 | -2998.16809 | 17.67 | 1.77 | 1.70 |
| r2SCAN-3c | -2997.38538 | -2997.36833 | 10.70 | 1.70 | 1.50 |

### Side-on quartet

| Method | RC (Hartree) | TS (Hartree) | ΔE (kcal/mol) | <S²> RC | <S²> TS |
|---|---:|---:|---:|---:|---:|
| HF/def2-SVP RIJCOSX def2/J | -2985.88043 | -2985.84826 | 20.19 | 6.62 | 6.98 |
| RI-MP2/def2-SVP | -2990.51794 | -2990.50015 | 11.16 | | |
| HF/def2-TZVP RIJCOSX def2/J | -2987.63067 | -2987.59641 | 21.50 | 6.66 | 7.00 |
| RI-MP2/def2-TZVP | -2993.34450 | -2993.32593 | 11.65 | | |
| HF/def2-QZVP RIJCOSX def2/J | -2987.73552 | -2987.70193 | 21.08 | 6.82 | 7.15 |
| RI-MP2/def2-QZVP | -2993.96842 | -2993.95088 | 11.01 | | |
| HF/cc-pVDZ RIJCOSX def2/J | -2987.30434 | -2987.27114 | 20.83 | 6.63 | 7.00 |
| RI-MP2/cc-pVDZ | -2992.07944 | -2992.05898 | 12.84 | | |
| HF/cc-pVTZ RIJCOSX def2/J | -2987.63521 | -2987.60142 | 21.21 | 6.79 | 7.13 |
| RI-MP2/cc-pVTZ | -2993.53580 | -2993.51709 | 11.74 | | |
| HF/cc-pVQZ RIJCOSX def2/J | -2987.72331 | -2987.68981 | 21.03 | 6.84 | 7.15 |
| RI-MP2/cc-pVQZ | -2994.04995 | -2994.03217 | 11.16 | | |
| B3LYP/def2-SVP RIJCOSX def2/J | -2995.47087 | -2995.44205 | 18.09 | 3.79 | 3.82 |
| B3LYP/def2-TZVP RIJCOSX def2/J | -2997.22146 | -2997.19150 | 18.80 | 3.78 | 3.83 |
| B3LYP/def2-QZVP RIJCOSX def2/J | -2997.34952 | -2997.32001 | 18.52 | 3.79 | 3.83 |
| B3LYP/cc-pVDZ RIJCOSX def2/J | -2996.84253 | -2996.81385 | 18.00 | 3.79 | 3.83 |
| B3LYP/cc-pVTZ RIJCOSX def2/J | -2997.22454 | -2997.19495 | 18.57 | 3.79 | 3.82 |
| B3LYP/cc-pVQZ RIJCOSX def2/J | -2997.32710 | -2997.29775 | 18.42 | 3.79 | 3.82 |
| B3LYP-D4/def2-SVP RIJCOSX def2/J | -2995.63475 | -2995.61252 | 13.95 | 3.78 | 3.83 |
| B3LYP-D4/def2-TZVP RIJCOSX def2/J | -2997.38534 | -2997.36197 | 14.67 | 3.79 | 3.83 |
| B3LYP-D4/def2-QZVP RIJCOSX def2/J | -2997.51340 | -2997.49048 | 14.38 | 3.79 | 3.82 |
| B3LYP-D4/cc-pVDZ RIJCOSX def2/J | -2997.00641 | -2996.98432 | 13.86 | 3.79 | 3.83 |
| B3LYP-D4/cc-pVTZ RIJCOSX def2/J | -2997.38842 | -2997.36542 | 14.43 | 3.79 | 3.82 |
| B3LYP-D4/cc-pVQZ RIJCOSX def2/J | -2997.49098 | -2997.46822 | 14.28 | 3.79 | 3.82 |
| r2SCAN-3c | -2997.38497 | -2997.36892 | 10.07 | 3.78 | 3.83 |

### Face-on doublet

| Method | RC (Hartree) | TS (Hartree) | ΔE (kcal/mol) | <S²> RC | <S²> TS |
|---|---:|---:|---:|---:|---:|
| HF/def2-SVP RIJCOSX def2/J | -2985.89660 | -2985.88100 | 9.79 | 5.73 | 6.18 |
| RI-MP2/def2-SVP | -2990.51145 | -2990.48524 | 16.45 | | |
| HF/def2-TZVP RIJCOSX def2/J | -2987.64155 | -2987.61960 | 13.77 | 5.49 | 6.19 |
| RI-MP2/def2-TZVP | -2993.35094 | -2993.30650 | 27.89 | | |
| HF/def2-QZVP RIJCOSX def2/J | -2987.74307 | -2987.72002 | 14.46 | 5.44 | 6.19 |
| RI-MP2/def2-QZVP | -2993.97564 | -2993.93765 | 23.84 | | |
| HF/cc-pVDZ RIJCOSX def2/J | -2987.31652 | -2987.29859 | 11.25 | 5.60 | 6.15 |
| RI-MP2/cc-pVDZ | -2992.07968 | -2992.04484 | 21.87 | | |
| HF/cc-pVTZ RIJCOSX def2/J | -2987.64375 | -2987.62147 | 13.98 | 5.43 | 6.18 |
| RI-MP2/cc-pVTZ | -2993.54304 | -2993.50519 | 23.75 | | |
| HF/cc-pVQZ RIJCOSX def2/J | -2987.73069 | -2987.70766 | 14.45 | 5.44 | 6.19 |
| RI-MP2/cc-pVQZ | -2994.05604 | -2994.01943 | 22.98 | | |
| B3LYP/def2-SVP RIJCOSX def2/J | -2995.46971 | -2995.43821 | 19.77 | 1.76 | 1.57 |
| B3LYP/def2-TZVP RIJCOSX def2/J | -2997.21925 | -2997.18681 | 20.36 | 1.76 | 1.56 |
| B3LYP/def2-QZVP RIJCOSX def2/J | -2997.34720 | -2997.31520 | 20.08 | 1.76 | 1.57 |
| B3LYP/cc-pVDZ RIJCOSX def2/J | -2996.84100 | -2996.81064 | 19.05 | 1.76 | 1.55 |
| B3LYP/cc-pVTZ RIJCOSX def2/J | -2997.22262 | -2997.19054 | 20.13 | 1.76 | 1.56 |
| B3LYP/cc-pVQZ RIJCOSX def2/J | -2997.32495 | -2997.29300 | 20.05 | 1.76 | 1.57 |
| B3LYP-D4/def2-SVP RIJCOSX def2/J | -2995.64311 | -2995.62049 | 14.20 | 1.76 | 1.57 |
| B3LYP-D4/def2-TZVP RIJCOSX def2/J | -2997.39266 | -2997.36909 | 14.79 | 1.76 | 1.56 |
| B3LYP-D4/def2-QZVP RIJCOSX def2/J | -2997.52064 | -2997.49747 | 14.54 | 1.76 | 1.57 |
| B3LYP-D4/cc-pVTZ RIJCOSX def2/J | -2997.39603 | -2997.37282 | 14.57 | 1.76 | 1.56 |
| B3LYP-D4/cc-pVQZ RIJCOSX def2/J | -2997.49828 | -2997.47525 | 14.45 | 1.76 | 1.57 |
| r2SCAN-3c | -2997.39211 | -2997.37396 | 11.39 | 1.68 | 1.36 |

### Face-on quartet

| Method | RC (Hartree) | TS (Hartree) | ΔE (kcal/mol) | <S²> RC | <S²> TS |
|---|---:|---:|---:|---:|---:|
| HF/def2-SVP RIJCOSX def2/J | -2985.87531 | -2985.83776 | 23.56 | 6.54 | 7.06 |
| RI-MP2/def2-SVP | -2990.52912 | -2990.49689 | 20.22 | | |
| HF/def2-TZVP RIJCOSX def2/J | -2987.62375 | -2987.58436 | 24.72 | 6.59 | 7.08 |
| RI-MP2/def2-TZVP | -2993.35534 | -2993.32403 | 19.65 | | |
| HF/def2-QZVP RIJCOSX def2/J | -2987.72820 | -2987.68989 | 24.04 | 6.76 | 7.22 |
| RI-MP2/def2-QZVP | -2993.97913 | -2993.94888 | 18.98 | | |
| HF/cc-pVDZ RIJCOSX def2/J | -2987.29707 | -2987.26067 | 22.84 | 6.54 | 7.07 |
| RI-MP2/cc-pVDZ | -2992.08712 | -2992.05583 | 19.63 | | |
| HF/cc-pVTZ RIJCOSX def2/J | -2987.62632 | -2987.58983 | 22.90 | 6.71 | 7.19 |
| RI-MP2/cc-pVTZ | -2993.54438 | -2993.51566 | 18.02 | | |
| HF/cc-pVQZ RIJCOSX def2/J | -2987.71405 | -2987.67797 | 22.64 | 6.77 | 7.24 |
| RI-MP2/cc-pVQZ | -2994.05796 | -2994.02976 | 17.69 | | |
| B3LYP/def2-SVP RIJCOSX def2/J | -2995.46931 | -2995.43762 | 19.89 | 3.78 | 3.82 |
| B3LYP/def2-TZVP RIJCOSX def2/J | -2997.21890 | -2997.18558 | 20.91 | 3.79 | 3.81 |
| B3LYP/def2-QZVP RIJCOSX def2/J | -2997.34669 | -2997.31386 | 20.60 | 3.79 | 3.81 |
| B3LYP/cc-pVDZ RIJCOSX def2/J | -2996.84043 | -2996.80887 | 19.80 | 3.79 | 3.81 |
| B3LYP/cc-pVTZ RIJCOSX def2/J | -2997.22211 | -2997.18916 | 20.68 | 3.79 | 3.81 |
| B3LYP/cc-pVQZ RIJCOSX def2/J | -2997.32438 | -2997.29164 | 20.54 | 3.79 | 3.81 |
| B3LYP-D4/def2-SVP RIJCOSX def2/J | -2995.64278 | -2995.61780 | 15.68 | 3.78 | 3.82 |
| B3LYP-D4/def2-TZVP RIJCOSX def2/J | -2997.39237 | -2997.36577 | 16.69 | 3.79 | 3.81 |
| B3LYP-D4/def2-QZVP RIJCOSX def2/J | -2997.52016 | -2997.49405 | 16.38 | 3.79 | 3.81 |
| B3LYP-D4/cc-pVDZ RIJCOSX def2/J | -2997.01390 | -2996.98906 | 15.59 | 3.79 | 3.81 |
| B3LYP-D4/cc-pVTZ RIJCOSX def2/J | -2997.39558 | -2997.36935 | 16.46 | 3.79 | 3.81 |
| B3LYP-D4/cc-pVQZ RIJCOSX def2/J | -2997.49787 | -2997.47182 | 16.35 | 3.79 | 3.81 |
| r2SCAN-3c | -2997.39089 | -2997.36975 | 13.27 | 3.79 | 3.81 |

## 3.2 DMRG

`results/dmrg/`

Classical reference energies for the Hamiltonians in
[Section 1.3](#13-active-space-hamiltonians), computed with an in-house
distributed multi-GPU DMRG implementation. These are the values against which
the quantum workflow is benchmarked. For CAS(11,10) an exact UCASCI reference is
also available and takes precedence.

One raw solver log per state and active space, named `{rc,ts}_{n}_{m}.out`. Each
log records the run configuration, a per-sweep history with energy and discarded
weight, and a final summary containing the bond-dimension ladder, the
extrapolation, ⟨S²⟩, natural-orbital occupations and an RDM energy consistency
check. The ladder is the part to parse:

```
Bond-dimension ladder:
      BD              E (Ha)          DW
    3000   -2987.4871247888   1.833e-04
    2500   -2987.4866720326   2.245e-04
    ...
```

### Best variational results

Energies at the largest bond dimension `D_max` of each ladder. ⟨S²⟩ is the
physical total-spin expectation value, evaluated with the alpha/beta MO overlap
(spin-pure doublet 0.75, broken-symmetry reference ≈ 1.75).

| Active space | `D_max` | E(RC) / Ha | ⟨S²⟩ RC | E(TS) / Ha | ⟨S²⟩ TS | ΔE‡ / kcal mol⁻¹ |
|---|---:|---:|---:|---:|---:|---:|
| CAS(11,10) | 200 | −2987.35539 | 3.09 | −2987.33167 | 3.65 | 14.9 |
| CAS(21,20) | 1200 | −2987.41960 | 3.27 | −2987.39349 | 4.19 | 16.4 |
| CAS(31,30) | 3000 | −2987.48712 | 4.29 | −2987.46514 | 4.63 | 13.8 |
| CAS(41,40) | 5000 | −2987.54391 | 4.53 | −2987.52557 | 5.15 | 11.5 |
| CAS(51,50) | 7000 | −2987.59998 | 4.76 | −2987.58133 | 5.65 | 11.7 |

### Extrapolated results

Linear fits of the energy against the discarded weight over the asymptotic part
of each ladder. `ε` is the maximum absolute residual of the fit; the barrier
uncertainty combines the RC and TS residuals in quadrature.

| Active space | E(δ→0) RC / Ha | ε / mHa | E(δ→0) TS / Ha | ε / mHa | ΔE‡ / kcal mol⁻¹ |
|---|---:|---:|---:|---:|---:|
| CAS(11,10) | −2987.35545 | 0.04 | −2987.33167 | <0.01 | 14.92 ± 0.02 |
| CAS(21,20) | −2987.42020 | 0.07 | −2987.39390 | 0.04 | 16.50 ± 0.05 |
| CAS(31,30) | −2987.49011 | 0.16 | −2987.46802 | 0.09 | 13.86 ± 0.12 |
| CAS(41,40) | −2987.55405 | 0.50 | −2987.53397 | 0.62 | 12.60 ± 0.50 |
| CAS(51,50) | −2987.61404 | 0.46 | −2987.59388 | 0.46 | 12.65 ± 0.41 |

The barrier first rises from CAS(11,10) to CAS(21,20), then falls and settles
near 12.6 kcal/mol. Small active spaces are clearly insufficient: between 20 and
50 orbitals the barrier moves by almost 4 kcal/mol, far beyond chemical
accuracy. The difference between the two largest active spaces is smaller than
the extrapolation uncertainty, which indicates the onset of convergence.

> **The extrapolation printed inside each `.out` file is not the value quoted
> above.** The solver fits every ladder point with `DW > 0`, whereas the
> published numbers exclude the lowest bond dimension of each ladder, where the
> discarded weight has already saturated while the energy is still moving. That
> point would bias the fit; dropping it keeps every residual below 0.7 mHa. To
> reproduce the table, parse the ladder and refit without its lowest `BD` row.

Two further properties matter when comparing against these numbers.

**Total energies include `E_core`.** RC and TS have different core energies, so
only total energies are subtractable. The barrier is `E(TS) − E(RC)`.

**⟨S²⟩ lies far above the spin-pure value.** The calculations constrain the spin
projection but not the total spin, and build on spin-contaminated unrestricted
orbitals, so the variational ground state acquires sizeable admixtures of higher
spin states, consistently more so for the TS than for the RC. A spin-adapted
solver will not reproduce these energies, and the difference is physical rather
than a convergence artifact.
