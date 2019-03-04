#!/usr/bin/env python

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2018 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
# This file is derived from software bearing the copyright notice in merlin4_copyright.txt

# Test that python, numpy and scipy are installed
from __future__ import print_function

import sys
print("Python version", sys.version)

try:
	import numpy
	print("numpy version", numpy.version.version)
except ImportError:
	print("ERROR: numpy not installed")
	exit(1)


try:
	import scipy
	import scipy.stats
	print("scipy version", scipy.version.version)
except ImportError:
	print("ERROR: scipy not installed")
	exit(1)




