import numpy as np
import psi4
import qm7_jk as jk
import time


n_trials = int(input("Enter number of trials: "))
basis_set = input("Enter basis set: ")

# Make sure we get the same random array
np.random.seed(0)
mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 104
""")

# Build a ERI tensor
basis = psi4.core.BasisSet.build(mol, target=basis_set)
mints = psi4.core.MintsHelper(basis)
I = np.array(mints.ao_eri())

# Symmetric random density
nbf = I.shape[0]
D = np.random.rand(nbf, nbf)
D = (D + D.T) / 2

# Einsum -------------------------------------------------------
print('-' * 50)
t_einsum = []
for i in range(n_trials):
    start = time.time()
    # Reference
    J_ref = np.einsum("pqrs,rs->pq", I, D, optimize=True)
    K_ref = np.einsum("prqs,rs->pq", I, D, optimize=True)
    end = time.time()
    t_einsum.append((end - start) * 100)
t_einsum_avg = sum(t_einsum) / len(t_einsum)
print('Einsum avg (%i runs): %.5f' % (n_trials, t_einsum_avg))

# Initialize density
J = np.random.rand(nbf, nbf)
K = np.random.rand(nbf,nbf)
D = D - np.diag(np.diag(D) * 0.5)

# QM7 JK -------------------------------------------------------
t_qm7 = []
for i in range(n_trials):
    start = time.time()
    jk.jk_numpy(I, D, J, K)
    end = time.time()
    t_qm7.append((end - start) * 100)
t_qm7_avg = sum(t_qm7) / len(t_qm7)
print('QM7 avg   (%i runs): %.5f' % (n_trials, t_qm7_avg))
print('-' * 50)
print('Ratio : ', t_qm7_avg / t_einsum_avg)
print('-' * 50)

# Make sure your implementation is correct
print("J is correct: %s" % np.allclose(J, J_ref))
print("K is correct: %s" % np.allclose(K, K_ref))
