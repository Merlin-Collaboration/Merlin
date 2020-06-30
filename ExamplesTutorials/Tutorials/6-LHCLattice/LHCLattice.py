#!/usr/bin/env python3

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2018 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
# This file is derived from software bearing the copyright notice in merlin4_copyright.txt

from __future__ import division, print_function
import numpy
from matplotlib import pyplot, gridspec
import os
import sys

data_filename = "tutorial6.out"
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
	
	grid = gridspec.GridSpec(2,1,height_ratios=[3,1])
	plt1 = pyplot.subplot(grid[0])
	plt1.set_ylabel("Beta [m]")

	beta_x, = plt1.plot(latticefunctions[:,0],latticefunctions[:,7], color='b', label='Beta x')
	beta_y, = plt1.plot(latticefunctions[:,0],latticefunctions[:,9], color='r', label='Beta y')
	plt1.legend(handles=[beta_x, beta_y],loc=1)

	plt2 = pyplot.subplot(grid[1])
	plt2.set_ylim([0,3.5])
	plt2.set_ylabel("Disp x [m]")	
	plt2.set_xlabel("Location [m]")
	disp_x, = plt2.plot(latticefunctions[:,0],(latticefunctions[:,11]/latticefunctions[:,15]), color='g', label='Disp x')
	plt2.legend(handles=[disp_x],loc=1)
		
	pyplot.show()
	
if all_good:
	print("Plotted lattice functions")
	os.remove(data_filename)
else:
	print("Failed to plot lattice functions. Could not read output file.")
	exit(1)

