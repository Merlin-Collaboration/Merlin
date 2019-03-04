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

data_filename = "tutorial4.out"
Particle1 = []
Particle2 = []
count = 1

data_paths = ["./build/","../build/", "../../build/","./build/output/","../build/output/","../../build/output/"]
all_good = False;
for data_path in data_paths:
	try:
		fname = os.path.join(data_path, data_filename)
		with open(fname, 'r') as file:
			for line in file :
				if count % 2 != 0 :
					Particle1.append(line.split())
				else :
					Particle2.append(line.split())
				count+=1
		print("Read test data:", fname)
	except IOError:
		continue#

	Particle1array = numpy.asarray(Particle1)
	Particle2array = numpy.asarray(Particle2)

	all_good = True;

	pyplot.xlabel("x [m]")
	pyplot.ylabel("px [m]")

	particle1, = pyplot.plot(Particle1array[:,2],Particle1array[:,3], color='b', label='1mm offset')
	particle2, = pyplot.plot(Particle2array[:,2],Particle2array[:,3], color='r', label='2mm offset')
	pyplot.legend(handles=[particle1, particle2],loc=1)

	pyplot.tight_layout()

	pyplot.show()
	
if all_good:
	print("Plotted lattice functions")
	os.remove(data_filename)
else:
	print("Failed to plot lattice functions. Could not read output file.")
	exit(1)

