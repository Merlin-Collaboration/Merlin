#!/usr/bin/env python3

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2018 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
# This file is derived from software bearing the copyright notice in merlin4_copyright.txt

from __future__ import division, print_function
import numpy
from matplotlib import pyplot
import os
import sys

data_filename = "tutorial5.out"
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

data_paths = ["../../","../../../","../../../../","../../../../../"]
all_good = False;
for data_path in data_paths:
	try:
		fname = os.path.join(data_path, data_filename)
		BPMcoord = numpy.loadtxt(fname)
		print("Read test data:", fname)
	except IOError:
		continue#

	all_good = True;

	pyplot.xlabel("x [m]")
	pyplot.ylabel("y [m]")

	centroid, = pyplot.plot(BPMcoord[:,0],BPMcoord[:,2], color='b', label='Bunch centroid')
	pyplot.legend(handles=[centroid],loc=1)
	pyplot.tight_layout()

	pyplot.show()
	
if all_good:
	print("Plotted lattice functions")
	os.remove(data_filename)
else:
	print("Failed to plot lattice functions. Could not read output file.")
	exit(1)

