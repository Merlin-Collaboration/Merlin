#!/usr/bin/env python
from __future__ import division
#from math import *
import numpy
#import scipy
#from matplotlib import pyplot
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


