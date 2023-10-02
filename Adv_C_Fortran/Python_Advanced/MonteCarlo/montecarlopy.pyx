# Approach 3: Cython
# Save the following code in a file named monte_carlo_cython.pyx
# Compile it using: cythonize -i monte_carlo_cython.pyx
# Import the resulting module
# monte_carlo_cython.pyx
import numpy as np
cimport numpy as np
cimport cython
 
def monte_carlo_cython(int nTrials):
    cdef int i, inside_count
    cdef double x, y
    inside_count = 0
    for i in range(nTrials):
        x = np.random.rand()
        y = np.random.rand()
        if x**2 + y**2 <= 1:
            inside_count += 1
    return 4 * inside_count / nTrials