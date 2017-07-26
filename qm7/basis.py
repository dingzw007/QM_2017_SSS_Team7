"""
Building the library of the basis functions.
"""
import psi4


def get_basis(mol, basis_set):
    """
    Builds mints for a molecule and its basis sets
    """
    bas = psi4.core.BasisSet.build(mol, target=basis_set)
    bas.print_out()

    mints = psi4.core.MintsHelper(bas)
    nbf = mints.nbf()

    if (nbf > 100):
        raise Exception("MOre than 100 basis functions")

    return mints
