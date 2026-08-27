"""Converge the UHF reference state of the RC (side-on doublet).

The doublet ground state of Compound I is a broken-symmetry solution: a triplet
iron-oxo unit coupled antiferromagnetically to a ligand radical. An SCF started
from a default guess is not guaranteed to find it. The published checkpoint
provides the starting orbitals, and the same checkpoint is used for both the RC
and the TS so that the two structures converge to the same electronic solution
and their energy difference is a meaningful reaction barrier.
"""

from pyscf import gto, scf

GEOMETRY = "../geometries/RC_side_on_doublet.xyz"
INIT_CHK = "RC_side_on_doublet_uhf.chk"

mol = gto.M(
    atom=GEOMETRY,
    basis="cc-pvdz",
    charge=0,
    spin=1,           # 2S = n_alpha - n_beta, i.e. a doublet
    verbose=4,
)

mf = scf.UHF(mol)
mf.conv_tol = 1e-8
mf.max_cycle = 200
mf.chkfile = "rc_uhf.chk"          # output; the published file is not overwritten

dm0 = mf.from_chk(INIT_CHK)              # initial guess from the published orbitals
mf = mf.newton()                         # second-order (Newton) UHF
mf.kernel(dm0=dm0)

assert mf.converged, "UHF did not converge"
print(f"UHF energy  = {mf.e_tot:.10f} Ha")
print(f"<S^2>, 2S+1 = {mf.spin_square()}")

# mf.mo_coeff, mf.mo_energy and mf.mo_occ now hold the reference orbitals used
# by every correlated calculation in this project.
