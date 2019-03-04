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

if "plot" in sys.argv:
	from matplotlib import pyplot

ref_run = False
if "ref_run" in sys.argv[3:]:
	ref_run = True


sim_filename = "basic_hollow_electron_lens_test_out.dat"
data_filename = "basic_hollow_electron_lens_test_out.dat"
dt = numpy.dtype({"names":["case","id","x","xp","y","yp"], "formats":["i8","f8","f8","f8","f8","f8"]})
tol = 1e-16

try:
	os.remove(sim_filename)
except OSError:
	pass
print("Running basic_hollow_electron_lens_test")
exe_name = os.path.join(os.path.dirname(__file__), "basic_hollow_electron_lens_test")

print("running:", [exe_name])
ret = subprocess.call([exe_name])

if ret != 0:
	print("Failed to run executable %s, error %s"%(exe_name, ret))
	exit(ret)



try:
	run_data = numpy.loadtxt(sim_filename, skiprows=1, dtype=dt)
except IOError:
	print("Failed to read", sim_filename)
	exit(1)
print("Read simulated data")

if not ref_run:
	# Find test data
	data_paths = ["../data", "data", "MerlinTests/data"]

	probs_raw = None
	for data_path in data_paths:
		try:
			fname = os.path.join(data_path, data_filename)
			ref_data = numpy.loadtxt(fname, skiprows=1, dtype=dt)
			header = open(fname).readline()
			print("Read test data:", fname)
		except IOError:
			continue

	if ref_data is None:
		print("Failed to open dist file:", data_filename)
		exit(1)



if "plot" in sys.argv:
	print(run_data.shape)
	for testcase in numpy.unique(run_data["case"]):
		run_data_t = run_data[run_data["case"] == testcase]

		fig, axs = pyplot.subplots(nrows=2, sharex=True)
		axs[0].plot(run_data_t["x"],run_data_t["xp"],'rx')
		axs[0].set_xlabel("x")
		axs[0].set_ylabel("xp")
		axs[0].set_title("Testcase %d"%testcase)
		axs[1].plot(run_data_t["y"],run_data_t["yp"],'rx')
		axs[1].set_xlabel("y")
		axs[1].set_ylabel("yp")

		if not ref_run:
			ref_data_t = ref_data[ref_data["case"] == testcase]
			axs[0].plot(ref_data_t["x"],ref_data_t["xp"],'k+')
			axs[1].plot(ref_data_t["y"],ref_data_t["yp"],'k+')
	pyplot.show()

if not ref_run:
	all_good = True

	if run_data.size != ref_data.size:
		all_good = False
		print("Run data size (%s) does not match reference data size (%s)"%(run_data.size, ref_data.size))

	for testcase in numpy.unique(ref_data["case"]):
		ref_data_t = ref_data[ref_data["case"] == testcase]
		run_data_t = run_data[run_data["case"] == testcase]

		for key in ["x","xp","y","yp"]:
			diff = numpy.fabs((ref_data_t[key] - run_data_t[key]))
			if diff.max() > tol:
				all_good = False
				worst_id = diff.argmax()
				print("Difference with respect to reference data case: %s, key: %s"%(testcase, key))
				print("Max diff: id=%s delta=%s ( %s vs %s)"%( worst_id, diff[worst_id], ref_data_t[key][worst_id], run_data_t[key][worst_id]))

	if all_good:
		print("Pass")
		os.remove(sim_filename)
	else:
		print("Fail")
		exit(1)
