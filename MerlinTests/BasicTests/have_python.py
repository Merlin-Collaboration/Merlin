#!/usr/bin/env python
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




