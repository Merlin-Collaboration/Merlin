#!/usr/bin/env python

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2018 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
# This file is derived from software bearing the copyright notice in merlin4_copyright.txt

from __future__ import division, print_function
import numpy
import subprocess
import os
from scipy.stats import chi2
import sys

# This test should fail occasionally (less than 0.1% of runs).
# cu50_test.cpp runs nparticles though a 50 cm graphite collimator, and histograms the x,xp,y,yp,dp distributions
# This is based on the tests in "TOOLS FOR PREDICTING CLEANING EFFICIENCY IN THE LHC" by R. Assmann et al, PAC2003
# cu50_test.py compares this to the probabilities in cu50_merlin_hist.dat which is generated previously by cu50_test.cpp
# Tests that we don't see a significant (p=1e-3) deviation from the reference (+/- 1%)
# tol, pval and nparticles choosen for fast run time, and low false negatives
# add the argument "plot" for a plot 
# add argument "ref_run" to generate a reference file

tol = 1e-2
pval = 1e-3

seed = 0

if len(sys.argv) > 1:
	seed = int(sys.argv[1]) # zero for seed from clock

npart = 10000
if len(sys.argv) > 2:
	npart = int(sys.argv[2])

scatter_mode_sixtrack = False
data_filename = "cu50_test_hist_m_10000000000.dat"
sim_filename = "cu50_test_hist_m_%s.dat"%npart
if "sixtrack" in sys.argv[3:]:
	scatter_mode_sixtrack = True
	data_filename = "cu50_test_hist_st_10000000000.dat"
	sim_filename = "cu50_test_hist_st_%s.dat"%npart

ref_run = False
if "ref_run" in sys.argv[3:]:
	ref_run = True

if not ref_run:
	# Find test data
	data_paths = ["../data", "data", "MerlinTests/data"]
	dt = numpy.dtype({"names":["bin","x","xp","y","yp","dp"], "formats":["i8","f8","f8","f8","f8","f8"]})
	probs_raw = None
	for data_path in data_paths:
		try:
			fname = os.path.join(data_path, data_filename)
			probs_raw = numpy.loadtxt(fname, skiprows=2, dtype=dt)
			header = open(fname).readline()
			print("Read test data:", fname)
		except IOError:
			continue

	if probs_raw is None:
		print("Failed to open dist file:", data_filename)
		exit(1)

	nparts_ref = 0
	for _s in header.split():
		if _s.startswith("npart"): nparts_ref = int(_s.partition("=")[2])

	probs = probs_raw
	for key in ["x","xp","y","yp","dp"]: probs[key] /= nparts_ref


	if "plot" in sys.argv:
		from matplotlib import pyplot
		fig, axes = pyplot.subplots(nrows=5, sharex=True)
		fig.set_size_inches(8,11)

for run_n in range(1):
	# multiple runs useful for visualising if differences are statistical or systematic
	
	# Run the test script
	try:
		os.remove(sim_filename)
	except OSError:
		pass
	print("Running cu50_test")
	exe_name = os.path.join(os.path.dirname(__file__), "cu50_test")
	args = [str(seed), str(npart)]
	if scatter_mode_sixtrack: args.append("sixtrack")
	print("running:", [exe_name]+args)
	ret = subprocess.call([exe_name]+args)

	if ret != 0:
		print("Failed to run executable %s, error %s"%(exe_name, ret))
		exit(ret)
	
	if ref_run:
		print("Made reference file:", sim_filename)
		exit(0)

	try:
		results_data = numpy.loadtxt(sim_filename, skiprows=2, dtype=dt)
		header = open(sim_filename).readline()
	except IOError:
		print("Failed to read", sim_filename)
		exit(1)
	print("Read simulated data")
	npart_a = 0
	for _s in header.split():
		if _s.startswith("npart"): npart_a = int(_s.partition("=")[2])
	print(npart_a, "particles")
	assert(npart == npart_a) # check right number of particles tracked

	expected = probs.copy()
	for key in ["x","xp","y","yp","dp"]: expected[key] *= npart

	all_good = True
	for key in ["x","xp","y","yp","dp"]:
		print("Testing",key)
		#mask by threshold
		bin_thresh = 50
		useful_bins = (expected[key] > bin_thresh)
		df = useful_bins.sum() - 1
		print("  Useful bins", useful_bins.sum())

		if df < 5:
			all_good == False
			print("  Too few particles tracked")
			continue

		with numpy.errstate(divide='ignore', invalid='ignore'): # hide warning divide by zero
			#sq_diff = ((results_data[key]-expected[key])/numpy.sqrt(expected[key]))**2
			sq_diff = ((results_data[key]-expected[key])/(numpy.sqrt(expected[key])+expected[key]*tol))**2
		sq_diff[expected[key] == 0] = 0 # remove nans, from when expected == zero
		chi_sq = (sq_diff * useful_bins).sum()
		
		print("  Chi^2 =", chi_sq)
		cdf = chi2.cdf(chi_sq, df)
		print("  cdf(chi^2)=", cdf)
		print("  percent of chi^2 higher even if true dist:", (1-cdf)*100)

		if (1-cdf) > pval:
			print(" ",key, "does not contain significant deviations")
			
		else:
			all_good = False
			print(" ",key, "does contain a significant deviation")
			worst_bins = numpy.argsort(sq_diff * useful_bins)
			print("  Worst 10 bins:")
			for i in worst_bins[-10:]:
				print("    Bin %5i: Expected %8.2f, got %8.2f (chi^2=%.2f)"%(i, expected[key][i], results_data[key][i], sq_diff[i]))

	if "plot" in sys.argv:
		for n, key in enumerate(["x","xp","y","yp","dp"]):
			axes[n].semilogy(expected["bin"], expected[key], "-k", label="expected")
			axes[n].semilogy(results_data["bin"], results_data[key], "-b", label="simulated")
			axes[n].annotate(key, xy=(0.15, 0.85), xycoords="axes fraction", horizontalalignment="center")
		pyplot.xlabel("Bin number")
		pyplot.xlim([0,expected["bin"][-1]])
		pyplot.tight_layout()

if "plot" in sys.argv:
	pyplot.show()

if all_good:
	print("Pass")
	os.remove(sim_filename)
else:
	print("Fail")
	exit(1)

