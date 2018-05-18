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
# lhc_collimation_test.cpp runs nparticles though 20 turns of the LHC ring
# lhc_collimation_test.py compares this to the probabilities in lhc_collimation_test_lossmap_10000000.dat which is generated previously by lhc_collimation_test.cpp
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
data_filename = "lhc_collimation_test_lossmap_10000000.dat"
sim_filename = "lhc_collimation_test_lossmap_%s.dat"%npart

ref_run = "ref_run" in sys.argv[3:]


if not ref_run:
	# Find test data
	data_paths = ["../data", "data", "MerlinTests/data"]
	dt = numpy.dtype({"names":["name","s","bin_start","loss","temp","len"], "formats":["a64","f8","f8","f8","f8","f8"]})
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
	for key in ["loss"]: probs[key] /= nparts_ref

	if "plot" in sys.argv:
		from matplotlib import pyplot
		fig, axes = pyplot.subplots(nrows=2, sharex=False)
		fig.set_size_inches(8,11)

def expand_lossmaps(lm1, lm2, add_1m_steps=False):
	"""lm1 and lm2 might not have the same primary keys
	
	create and return new arrays that have all primary keys, and fill data key with 0 where needed.
	"""
	assert(lm1.dtype.names == lm2.dtype.names)
	ic = lambda _a,_b:abs(_a-_b)<1e-4
	pk = "s" # primary key
	dk = "loss" # data key
	all_keys = set(lm1[pk]) | set(lm2[pk])
	if add_1m_steps:
		all_keys |= set(numpy.arange(0,26000,1))
	
	all_keys = sorted(list(all_keys))
	nlm1 = numpy.zeros([len(all_keys)], dtype=lm1.dtype)
	nlm2 = numpy.zeros([len(all_keys)], dtype=lm1.dtype)
	
	n1, n2 = 0, 0
	for n,k in enumerate(all_keys):
		nlm1[n][pk] = k
		nlm2[n][pk] = k
		
		if n1<len(lm1) and ic(lm1[n1][pk],k):
			nlm1[n] = lm1[n1]
			n1 += 1

		if n2<len(lm2) and ic(lm2[n2][pk],k):
			nlm2[n] = lm2[n2]
			n2 += 1
	return nlm1, nlm2

for run_n in range(1):
	# multiple runs useful for visualising if differences are statistical or systematic
	
	# Run the test script
	try:
		os.remove(sim_filename)
	except OSError:
		pass
	print("Running lhc_collimation_test")
	exe_name = os.path.join(os.path.dirname(__file__), "lhc_collimation_test")
	args = [str(seed), str(npart)]

	print("running:", [exe_name]+args)
	ret = subprocess.call([exe_name]+args)
	ret=0

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
	for key in ["loss"]: expected[key] *= npart

	expected, results_data = expand_lossmaps(expected, results_data)

	all_good = True
	for key in ["loss"]:
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
		expected, results_data = expand_lossmaps(expected, results_data, add_1m_steps=1)
		for n, key in enumerate(["loss"]):
			axes[2*n].semilogy(expected[key], "-k", label="expected")
			axes[2*n].semilogy(results_data[key], "-b", label="simulated")
			axes[2*n].annotate(key, xy=(0.15, 0.85), xycoords="axes fraction", horizontalalignment="center")
			
			mask = numpy.logical_and(expected["s"]>19500, expected["s"]< 20500)
			axes[2*n+1].semilogy(expected[key][mask], "-k", label="expected")
			axes[2*n+1].semilogy(results_data[key][mask], "-b", label="simulated")
			
		pyplot.xlabel("Bin number")
		pyplot.tight_layout()

if "plot" in sys.argv:
	pyplot.show()

if all_good:
	print("Pass")
	os.remove(sim_filename)
else:
	print("Fail")
	exit(1)

