from qm7 import mol
import psi4


def test_molecule_should_be_read_correctly():
    water = mol.molecule("tests/water")
    assert water.natom() == 3
    assert type(water) == psi4.core.Molecule
