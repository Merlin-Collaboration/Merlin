#!/usr/bin/env python

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2018 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
# This file is derived from software bearing the copyright notice in merlin4_copyright.txt

from __future__ import division, print_function
import numpy
from matplotlib import pyplot
import os
import sys

data_filenames = ["tutorial3a.out","tutorial3b.out"]
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

data_paths = ["../../","../../../","../../../../","../../../../../"]
all_good = False;
for data_path in data_paths:
	try:
		fname = os.path.join(data_path, data_filenames[0])
		latticefunctions1 = numpy.loadtxt(fname)
		fname2 = os.path.join(data_path, data_filenames[1])
		latticefunctions2 = numpy.loadtxt(fname2)
		print("Read test data:", fname)
		print("Read test data:", fname2)
	except IOError:
		continue#

	all_good = True;
	pyplot.subplot(2,1,1)
	pyplot.ylabel("Beta x [m]")

	beta_x1, = pyplot.plot(latticefunctions1[:,0],latticefunctions1[:,7], color='b', label='Untouched')
	beta_x2, = pyplot.plot(latticefunctions2[:,0],latticefunctions2[:,7], color='b', linestyle='--', label='Manipulated')
	pyplot.legend(handles=[beta_x1, beta_x2],loc=1)

	pyplot.subplot(2,1,2)
	pyplot.ylabel("Beta y [m]")
	pyplot.xlabel("Location [m]")
	beta_x1, = pyplot.plot(latticefunctions1[:,0],latticefunctions1[:,9], color='r', label='Untouched')
	beta_x2, = pyplot.plot(latticefunctions2[:,0],latticefunctions2[:,9], color='r', linestyle='--', label='Manipulated')
	pyplot.legend(handles=[beta_x1, beta_x2],loc=1)


	pyplot.show()
	
if all_good:
	print("Plotted lattice functions")
	os.remove(data_filenames)
else:
	print("Failed to plot lattice functions. Could not read output file.")
	exit(1)

