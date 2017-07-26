"""
Building the library of the molecules
"""
import psi4


def molecule(mol_file):
    """
    Reads a z matix file and builds the geometry 
    """
    with open(mol_file, "r") as mf:
        mol_lines = mf.readlines() 
    mol = psi4.geometry("".join(mol_lines))

    #Build a molecule
    mol.update_geometry()
    mol.print_out()

    return mol



