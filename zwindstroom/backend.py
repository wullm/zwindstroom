"""
Interface with the shared C library via ctypes.
"""

import ctypes
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
engine = ctypes.CDLL(dir_path + "/../build/libzwindstroom.so")
