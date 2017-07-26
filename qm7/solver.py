"""
This is to solve the converged charge matrix
"""
import psi4
import numpy as np


def solve(mints, nel, mol, n_iterations, e_conv, d_conv, damp_value, damp_start):
    """
    SCF
    """
    V = np.array(mints.ao_potential())
    T = np.array(mints.ao_kinetic())

    H = T + V
    S = np.array(mints.ao_overlap())
    g = np.array(mints.ao_eri())

    A = mints.ao_overlap()
    A.power(-0.5, 1.e-14)
    A = np.array(A)

    def diag(F, A):
        Fp = A.T @ F @ A
        eps, cp = np.linalg.eigh(Fp)
        C = A @ cp
        return eps, C
    eps, C = diag(H, A)
    Cocc = C[:, :nel]
    D = Cocc @ Cocc.T
    print(' # |      E_total     |    E_diff   |   Grad_rms  |')
    E_old = 0.0
    F_old = None
    # Update Fock matrix by calculating new density
    for i in range(n_iterations):
        # Fock matrix -> F_pq = H_pq + 2 * g_pqrs D_rs - g_pqrs D_rs
        J = np.einsum("pqrs,rs->pq", g, D)  # same as np.sum(g * D, axis=(2, 3))
        K = np.einsum("prqs,rs->pq", g, D)

        F_new = H + 2.0 * J - K   # Fock matrix

        # Damping
        if i > damp_start:
            F = (damp_value) * F_old + (1 - damp_value) * F_new
        else:
            F = F_new
        F_old = F_new

        # Build the AO gradient to determine convergence
        grad = F @ D @ S - S @ D @ F
        grad_rms = np.mean(grad ** 2) ** 0.5

        E_electric = np.sum((F + H) * D)
        E_total = E_electric + mol.nuclear_repulsion_energy()

        E_diff = E_total - E_old
        E_old = E_total
        print("%2i | % 16.12f | % 8.4e | % 8.4e |" % (i + 1, E_total, E_diff, grad_rms))

        if (E_diff <= e_conv) and (grad_rms <= d_conv):
            break

        eps, C = diag(F, A)
        Cocc = C[:, :nel]
        D = Cocc @ Cocc.T

    print("SCF has finished!\n")
