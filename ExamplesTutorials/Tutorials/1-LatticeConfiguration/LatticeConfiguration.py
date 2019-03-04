#!/usr/bin/env python

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2018 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
# This file is derived from software bearing the copyright notice in merlin4_copyright.txt

from __future__ import division, print_function
import numpy
from matplotlib import pyplot, gridspec
import os
import sys

data_filename = "tutorial1.out"
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

data_paths = ["../../","../../../","../../../../","../../../../../"]
all_good = False;
for data_path in data_paths:
	try:
		fname = os.path.join(data_path, data_filename)
		latticefunctions = numpy.loadtxt(fname)
		print("Read test data:", fname)
	except IOError:
			continue#

	all_good = True;
	
	pyplot.ylabel("Beta [m]")
	pyplot.xlabel("Location [m]")

	beta_x, = pyplot.plot(latticefunctions[:,0],latticefunctions[:,7], color='b', label='Beta x')
	pyplot.legend(handles=[beta_x],loc=1)
		
	pyplot.show()
	
if all_good:
	print("Plotted lattice functions")
	os.remove(data_filename)
else:
	print("Failed to plot lattice functions. Could not read output file.")
	exit(1)

