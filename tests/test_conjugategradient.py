import numpy as np
import qm7.conjugategradient as cg


def test_cg():
    A = np.array([[1, 2],
                 [2, 1]])
    b = np.array([3, 3])
    x_init = np.array([0.01, 0.01])
    error = cg.conjugate_gradient(A, b, x_init) - np.array([1, 1])
    assert np.linalg.norm(error) < 1e-8	
