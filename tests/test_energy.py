"""
this is to test the eneryg of H2O with sto-3g basis

"""

from qm7.scf import SCF
import psi4
import argparse
import numpy as np


def test_energy():
    args = argparse.Namespace()
    args.molecule = 'tests/water'
    args.basis = 'sto-3g'
    args.nocc = 5
    args.max_iter = 50
    args.econv = 1e-6
    args.dconv = 1e-6
    args.dampvalue = 0.2
    args.dampstart = 5
    args.damping = False
    args.DIIS = True

    scf = SCF(args)
    scf.initialize()
    scf.solve()

    psi4.set_options({"scf_type":"pk"})
    psi4_energy = psi4.energy("SCF/sto-3g", molecule=scf.mol)
    assert np.allclose(scf.Etotal, psi4_energy)
