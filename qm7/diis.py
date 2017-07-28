import numpy as np


def diis(F_list, DIIS_RES):
    """
    Return Fock matrix using Fock matrix list and DIIS residuals
    """
    B_dim = len(F_list) + 1
    B = np.empty((B_dim, B_dim))
    B[-1, :] = -1
    B[:, -1] = -1
    B[-1, -1] = 0
    for i in range(len(F_list)):
        for j in range(len(F_list)):
            B[i, j] = np.einsum('ij,ij->', DIIS_RES[i], DIIS_RES[j])

    # Build RHS of Pulay equation
    rhs = np.zeros((B_dim))
    rhs[-1] = -1

    # Solve Pulay equation for c_i's with NumPy
    coeff = np.linalg.solve(B, rhs)

    # Build DIIS Fock matrix
    F = np.zeros_like(F_list[-1])
    for x in range(coeff.shape[0] - 1):
        F += coeff[x] * F_list[x]

    return F
