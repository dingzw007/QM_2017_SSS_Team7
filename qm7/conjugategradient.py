"""
This is a generalized conjugate gradient solver

"""
import numpy as np

def conjugate_gradient(A, b, x_init):
    """
    This is a generalized conjugate gradient solver.
    You'll need to feed it the following:
    1) an Array A that is symmetric and positive definite.
    2) a Vector b 
    3) an initial guess x_init, which has the same dimensions as b.
    
    """
    tol = 1e-12 #Set tolerance to 1e-14, edit as per requirements
    iter_max = 1000 #Set maximum iterations to 1000, edit as per requirements

    #Init

    x = x_init
    g = b - np.dot(A, x)
    k = g.copy()

    #iterate

    for iteration in range(iter_max):


        p = np.vdot(g.T, g)
        q = np.dot(k.T, A)
        a = p/np.vdot(q, k)
        x = x + k * a
        h = g - np.dot(A * a, k)

        print("Residual %12.10e" % np.linalg.norm(h))
        if np.linalg.norm(h) < tol:
            return x
        b = np.vdot(h.T, h)/np.vdot(g.T, g)
        k = h + b * k
        g = h
    return x
    
    
