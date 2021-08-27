# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:03:14 2019

@author: kst

Setup op file, for "compiling" cython scripts. 
run in the terminal as:
python setup.py build_ext --inplace
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy 

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension('GarblingFreeXOR_c', 
                             ["GarblingFreeXOR_c.pyx"],
                             include_dirs = [numpy.get_include()])]
)
