"""Code to dress a Hamiltonian in FCIDUMP format with orbital rotations defined in a json file.

The json file should contain a list of tuples (j, i, theta) where j and i are the indices of the orbitals to be rotated and theta is the rotation angle.
The code reads the Hamiltonian from the FCIDUMP file, applies the orbital rotations, and saves the dressed Hamiltonian in a new FCIDUMP file.

Usage:
    fcidump_filename = "n2_6e6o_ccpvdz.FCIDUMP"
    rotation_filename = "mock_orbital_rotations.json"
    dress_fcidump(fcidump_filename, rotation_filename)
    
Output:
    The dressed Hamiltonian is saved in a new FCIDUMP file named "dressed_" + fcidump_filename.
"""

import json
import os
from itertools import combinations

import numpy as np
from pyscf import ao2mo, gto, lib, scf
from pyscf.tools import fcidump
from scipy.linalg import expm


def _oo_params(filename: str, nmo: int) -> list[tuple[int, int, list[float]]]:
    """Creation of the orbital rotation parameters from a json file.
    The json file should contain a list of tuples (j, i, theta) where j and i are the indices of the orbitals to be rotated and theta is the rotation angle.
    The function returns one list of tuples, each tuple containing the indices of the orbitals to be rotated and the rotation angle.
    Parameters
    ----------
    filename : str
        The name of the json file containing the orbital rotation parameters.
    nmo : int
        The number of spatial orbitals in the system.
    Returns
    -------
    list[tuple[int, int, list[float]]]
        A list of tuples, each tuple containing the indices of the orbitals to be rotated and the rotation angle.
    """
    # check if the file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File {filename} not found.")

    with open(filename, "r") as f_json:
        a = json.load(f_json)["act_rot"]
    oo_params = []
    for j, i, theta in a:
        if i < nmo and j < nmo:
            oo_params.append((j, i, theta))
        else:
            raise ValueError(f"Invalid rotation index: i={i}, j={j} for nom={nmo}")
    print(f"Loaded {len(oo_params)} orbital rotation parameters from {filename}")
    assert len(oo_params) == nmo * (nmo - 1) // 2, (
        "Number of rotations should be nmo choose 2"
    )
    return oo_params


def _define_U(gen_rot: list[tuple[int, int, float]], nmo: int) -> np.ndarray:
    """
    Define the unitary transformation U from the orbital rotation parameters.
    U is defined as the product of the exponentials of the generators of the orbital rotations.
    Parameters
    ----------
    gen_rot : list[tuple[int, int, float]]
        A list of tuples, each tuple containing the indices of the orbitals to be rotated and the rotation angle.
    nmo : int
        The number of spatial orbitals in the system.
    Returns
    -------
    np.ndarray
        The unitary transformation U as a numpy array.
    """

    T_R = []

    print("Generators:", len(gen_rot))
    for j, i, r in gen_rot:
        t = np.zeros((nmo, nmo), dtype=float)
        t[i][j] = -r
        t[j][i] = r

        T_R.append(t)

    U_rot = np.eye(nmo, dtype=float)
    for t in T_R:
        U_rot = expm(t) @ U_rot

    print("Unitary U:", np.allclose(U_rot @ U_rot.T, np.eye(nmo)))
    return U_rot


def _apply_U(
    h1e: np.ndarray,
    h2e: np.ndarray,
    U_rot: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply the unitary transformation to the Hamiltonian.
    """

    h1e_new = U_rot.T @ h1e @ U_rot

    h2e_new = lib.einsum("pqrs,pi,qj,rk,sl->ijkl", h2e, U_rot, U_rot, U_rot, U_rot)

    return h1e_new, h2e_new


def dress_fcidump(fcidump_filename: str, rotation_filename: str) -> None:
    """
    Dress the Hamiltonian in the FCIDUMP file with the orbital rotations defined in the rotation file.
    The dressed Hamiltonian is saved in a new FCIDUMP file.
    Parameters
    ----------
    fcidump_filename : str
        The name of the FCIDUMP file containing the Hamiltonian to be dressed.
    rotation_filename : str
        The name of the json file containing the orbital rotation parameters.
    Returns
    -------
    None
    """

    # Load the Hamiltonian from the FCIDUMP file
    data = fcidump.read("n2_6e6o_ccpvdz.fcidump")
    h1e = data["H1"]  # 1-electron integrals
    _h2e = data["H2"]  # 2-electron integrals
    ecore = data["ECORE"]  # Core nuclear repulsion energy
    nmo = data["NORB"]  # Number of molecular orbitals
    nelec = data["NELEC"]  # Number of electrons
    ms2 = data["MS2"]  # Spin multiplicity
    orbsym = data["ORBSYM"]  # Orbital symmetries

    # HERE YOU RESTORE BACK TO THE RIGHT 4-RANK tensor
    h2e = ao2mo.restore(1, _h2e, nmo)  # restore to full 4-index tensor
    print("\n------------FCIDUMP DATA------------")
    print(f"Loaded Hamiltonian from: {fcidump_filename}:")
    print(f"Number of molecular orbitals: {nmo}")
    print(f"Number of electrons: {nelec}")
    print(f"Core nuclear repulsion energy: {ecore}")
    print(f"Shape of 1-electron integrals: {h1e.shape}")
    print(f"Shape of 2-electron integrals: {h2e.shape}")

    print("\n------------DRESSED HAMILTONIAN------------")
    # Load the orbital rotation parameters from the json file
    oo_params = _oo_params(rotation_filename, nmo)

    # Define the unitary transformation U from the orbital rotation parameters
    U_rot = _define_U(oo_params, nmo)

    # Apply the unitary transformation to the Hamiltonian
    h1e_new, h2e_new = _apply_U(h1e, h2e, U_rot)

    assert h1e_new.shape == h1e.shape, (
        "Shape of 1-electron integrals has changed after dressing"
    )
    assert h2e_new.shape == h2e.shape, (
        "Shape of 2-electron integrals has changed after dressing"
    )

    print("h1e changed:", not np.allclose(h1e, h1e_new))
    print("h2e changed:", not np.allclose(h2e, h2e_new))
    # Save the dressed Hamiltonian in a new FCIDUMP file
    fcidump.from_integrals(
        filename="dressed_" + fcidump_filename,
        h1e=h1e_new,
        h2e=h2e_new,
        nuc=ecore,
        nmo=nmo,
        nelec=nelec,
        ms=ms2,
        orbsym=orbsym,
    )

    print(f"Dressed Hamiltonian saved in: dressed_{fcidump_filename}")


if __name__ == "__main__":
    from pyscf import mcscf

    mol = gto.M(
        atom="N 0 0 0; N 0 0 1.098",
        basis="cc-pvdz",
        verbose=2,
    )

    ncas = 6
    nelec = 6

    mf = scf.RHF(mol).run()
    mycas_ref = mcscf.CASCI(mf, ncas, nelec)

    # Save the Hamiltonian in FCIDUMP format
    fcidump.from_mcscf(mycas_ref, filename="n2_6e6o_ccpvdz.FCIDUMP")

    # mock up of the orbital rotation parameters
    rng = np.random.default_rng(seed=42)
    rotation_indexes = combinations(range(ncas), 2)
    active_rotations = []

    for i, j in rotation_indexes:
        active_rotations.append((j, i, rng.uniform(-1, 1)))
        print(f"Phi_[{j},{i}] = {active_rotations[-1][2]:.4f}")

    assert len(active_rotations) == ncas * (ncas - 1) // 2, (
        "Number of rotations should be ncas choose 2"
    )

    # mock up of the json file containing the orbital rotation parameters
    # these parameters are generally retrieved directly from the ADAPT-VMPE circuits
    rotation_filename = "mock_orbital_rotations.json"
    with open(rotation_filename, "w") as f_json:
        json.dump({"act_rot": active_rotations}, f_json, indent=4)

    # Dress the Hamiltonian in the FCIDUMP file with the orbital rotations defined in the rotation file
    dress_fcidump("n2_6e6o_ccpvdz.FCIDUMP", rotation_filename)
