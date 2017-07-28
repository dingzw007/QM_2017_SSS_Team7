"""
Self consistent field
"""
import psi4
import numpy as np
from qm7.basis import get_basis
from qm7.mol import molecule


class SCF:
    """ Self consistent field class """
    def __init__(self, cli_args):
        self.mol = molecule(cli_args.molecule)
        self.basis = get_basis(self.mol, cli_args.basis)
        self.nel = cli_args.nob
        self.n_iterations = cli_args.iteration
        self.e_conv = cli_args.econv
        self.d_conv = cli_args.dconv
        self.damp_value = cli_args.dampvalue
        self.damp_start = cli_args.dampstart
        self.Etotal = 0.0

    def initialize(self):
        """
        Initialize SCF
        """
        self.V = np.array(self.basis.ao_potential())
        self.T = np.array(self.basis.ao_kinetic())

        self.H = self.T + self.V
        self.S = np.array(self.basis.ao_overlap())
        self.g = np.array(self.basis.ao_eri())

        self.A = self.basis.ao_overlap()
        self.A.power(-0.5, 1.e-14)
        self.A = np.array(self.A)

        self.eps, self.C = self.diag(self.H, self.A)
        self.Cocc = self.C[:, :self.nel]
        self.D = self.Cocc @ self.Cocc.T

    def diag(self, F, A):
        Fp = A.T @ F @ A
        eps, cp = np.linalg.eigh(Fp)
        C = A @ cp
        return eps, C

    def solve(self):
        """
        Solve SCF
        """
        print(' # |      E_total     |    E_diff   |   Grad_rms  |')
        E_old = 0.0
        F_old = None
        for i in range(self.n_iterations):
            J = np.einsum("pqrs,rs->pq", self.g, self.D)
            K = np.einsum("prqs,rs->pq", self.g, self.D)

            F_new = self.H + 2.0 * J - K

            # Damping
            if i > self.damp_start:
                self.F = (self.damp_value) * F_old + (1 - self.damp_value) * F_new
            else:
                self.F = F_new
            F_old = F_new

            # Build the AO gradient to determine convergence
            grad = self.F @ self.D @ self.S - self.S @ self.D @ self.F
            grad_rms = np.mean(grad ** 2) ** 0.5

            E_electric = np.sum((self.F + self.H) * self.D)
            E_total = E_electric + self.mol.nuclear_repulsion_energy()

            E_diff = E_total - E_old
            E_old = E_total
            print("%2i | % 16.12f | % 8.4e | % 8.4e |" % (i + 1, E_total, E_diff, grad_rms))

            if (E_diff <= self.e_conv) and (grad_rms <= self.d_conv):
                self.Etotal = E_total
                break

            self.eps, self.C = self.diag(self.F, self.A)
            self.Cocc = self.C[:, :self.nel]
            self.D = self.Cocc @ self.Cocc.T

        print("SCF has finished!\n")
