#!/usr/bin/env python

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2018 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING

from __future__ import division, print_function
import os, sys
import numpy
from scipy import integrate
import scipy.stats

import pygsl.rng

nbins = 2000
outfile_path = "random_test_gendata.dat"

def find_data_file(fname):
	paths = ["MerlinTests/data/", "data/", "../data/"]

	for path in paths:
		if os.path.exists(os.path.join(path, fname)):
			return os.path.join(path, fname)
	print("Could not find data file '%s'. Try running from cmake build directory, or the executable directory"% fname)
	sys.exit(1)

random_setup_fname = find_data_file("random_test_setup.tfs")
dt = numpy.dtype({"names":["name", "dist", "var1", "var2", "var3", "range_min", "range_max"], "formats":["U32", "U32", "f", "f", "f", "f", "f"]})
random_setup = numpy.loadtxt(random_setup_fname, skiprows=2, dtype=dt)

print(random_setup)

test_hists = {}
for mn, mode in enumerate(random_setup):
	print(mode["name"])

	var1, var2, var3 = mode["var1"], mode["var2"], mode["var3"]

	hist = numpy.zeros([nbins+2])
	dist_type = "continuous"

	if mode["dist"] == "uniform":
		pdf  = scipy.stats.uniform(loc=var1, scale=var2-var1).pdf
	elif mode["dist"] == "normal":
		pdf  = scipy.stats.norm(loc=var1, scale=var2**0.5).pdf
	elif mode["dist"] == "normalc":
		a = var1 - (var2*var3)
		b = var1 + (var2*var3)
		pdf  = scipy.stats.truncnorm(a=a, b=b, loc=var1, scale=var2).pdf
	elif mode["dist"] == "poisson":
		dist_type = "discrete"
		pdf  = scipy.stats.poisson(mu=var1) # not actually the pdf, but we deal with that later
	elif mode["dist"] == "landau":
		pdf  = pygsl.rng.landau_pdf
	else:
		print("WARN: using wrong dist")
		pdf  = scipy.stats.uniform(loc=1, scale=2)

	if dist_type == "continuous":
		bin_min, bin_max = mode["range_min"],mode["range_max"]
		bin_width = (bin_max - bin_min) / (nbins)

		hist[0] = integrate.quad(pdf, -numpy.inf, bin_min)[0]
		hist[-1] = integrate.quad(pdf, bin_max, numpy.inf)[0]
		for n in range(nbins):
			this_bin_start = bin_min + n*bin_width
			this_bin_end = bin_min + (n+1)*bin_width

			hist[n+1] = integrate.quad(pdf, this_bin_start ,this_bin_end )[0]
	else:
		dist = pdf # deal with there not being a pdf for discrete
		bin_min = int(mode["range_min"])
		bin_max = int(mode["range_max"])
		hist -= 1
		#bin_max = bin_min + nbins
		hist[0] = dist.cdf(bin_min)
		hist[bin_max-bin_min+1] = 1 - dist.cdf(bin_max-1)

		hist[1:bin_max-bin_min+1] = dist.pmf(numpy.arange(bin_min, bin_max, dtype=int))

	print("Sum should be 1: ", hist[hist>=0].sum())
	test_hists[mode["name"]] = hist



with open(outfile_path, "w") as outfile:
	outfile.write("#generated with BasicTests/random_test_gendata.py\n")
	outfile.write(" ".join(random_setup['name'])+"\n")
	for n in range(nbins+2):
		print(n, file=outfile, end=" ")
		for dname in random_setup['name']:
			print (test_hists[dname][n], file=outfile, end=" ")
		print("", file=outfile)
