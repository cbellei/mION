import sys

try:
   import numpypy as np    # for compatibility with numpy in pypy
except:
   import numpy as np      # if using numpy in cpython

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    nf = len(b)     # number of equations
    a = np.append(0.,a)
    c = np.append(c,0.)
    for it in xrange(1, nf):
        m = a[it]/b[it-1]
        b[it] = b[it] - m*c[it-1] 
        d[it] = d[it] - m*d[it-1]
        	    
    x = np.zeros(len(d))
    x[-1] = d[-1]/b[-1]

    for il in xrange(nf-2, -1, -1):
        x[il] = (d[il]-c[il]*x[il+1])/b[il]

    return x