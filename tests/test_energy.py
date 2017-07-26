"""
this is to test the eneryg of H2O with sto-3g basis

"""

from qm7.mol import molecule
from qm7.basis import get_basis
from qm7.solver import solve
import psi4


def test_energy():
    mol = molecule("tests/water")
    basis = get_basis(mol,"sto-3g")
    test_energy = solve(basis,5,mol)
    psi4.set_options({"scf_type":"pk"})
    psi4_energy = psi4.energy("SCF/sto-3g",molecule=mol)
    print(test_energy)
    print(psi4_energy)
    assert  abs(test_energy-psi4_energy)<0.0001  
