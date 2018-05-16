#!/usr/bin/env python

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2018 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
# This file is derived from software bearing the copyright notice in merlin4_copyright.txt

from __future__ import division
import numpy
from scipy import integrate
import pygsl.rng


nbins = 2600
bin_min = -5
bin_max = 20
bin_width = (bin_max - bin_min) / (nbins)

hist = numpy.zeros([nbins+2])
# first and last bin catch outsiders
hist[0] = integrate.quad(pygsl.rng.landau_pdf, -numpy.inf, bin_min)[0]
hist[-1] = integrate.quad(pygsl.rng.landau_pdf, bin_max, numpy.inf)[0]


for n in xrange(nbins):
	this_bin_start = bin_min + n*bin_width
	this_bin_end = bin_min + (n+1)*bin_width

	hist[n+1] = integrate.quad(pygsl.rng.landau_pdf, this_bin_start ,this_bin_end )[0]


with open("../data/landau_test_dist.dat", "w") as outfile:
	outfile.write("#bin, norm_count\n")
	outfile.write("#generated with BasicTests/landau_test_gendata.py\n")
	for n in xrange(nbins+2):
		print >>outfile, n, hist[n]


