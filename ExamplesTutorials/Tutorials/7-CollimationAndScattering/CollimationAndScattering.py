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

data_filename = "tutorial7.out"
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

data_paths = ["", "../", "../../","../../../","../../../../","../../../../../"]
all_good = False;
for data_path in data_paths:
	try:
		fname = os.path.join(data_path, data_filename)
		data = numpy.loadtxt(fname,skiprows=1,usecols=(1,2,3))
		print("Read test data:", fname)
	except IOError:
		continue

	all_good = True;

	pyplot.xlabel("Location [m]")
	pyplot.ylabel("Loss Count")
	pyplot.yscale('log')
	
	pyplot.hist(data[:,0], weights=data[:,2], bins = 120)
	pyplot.legend(['Lost Particles'],loc=1)
	pyplot.tight_layout()

	pyplot.show()
	break

if all_good:
	print("Plotted lattice functions")
	os.remove(fname)
else:
	print("Failed to plot lattice functions. Could not read output file.")
	exit(1)

