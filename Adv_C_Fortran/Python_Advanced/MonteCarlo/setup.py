"""
 Run this code with the command:
    > python3 setup.py build_ext --inplace
"""
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize("montecarlopy.pyx"),
    include_dirs=[numpy.get_include()]
)   
