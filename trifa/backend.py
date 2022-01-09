"""
Contains methods to set up and execute the cosmological
and fluid calculations
"""

import ctypes
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
engine = ctypes.CDLL(dir_path + "/../build/lib3fa.so")
