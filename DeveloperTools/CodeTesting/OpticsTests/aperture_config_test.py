#!/usr/bin/env python3
from __future__ import division, print_function

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2018 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING

# Script to run and check output from aperture_config_test.cpp

import os
import sys
import subprocess
from math import *
import numpy
import itertools

if sys.argv.index("--exe"):
	i = sys.argv.index("--exe")
	sys.argv.pop(i)
	exe_name = sys.argv.pop(i)
else:
	exe_name = os.path.join(os.path.dirname(__file__), "aperture_config_test")

def find_data_file(fname):
	paths = ["DeveloperTools/CodeTesting/data/", "data/", "../data/"]

	for path in paths:
		if os.path.exists(os.path.join(path, fname)):
			return os.path.join(path, fname)
	print("Could not find data file '%s'. Try running from cmake build directory, or the executable directory"% fname)
	sys.exit(1)

def main():
	print("Running aperture_config_test")
	args = []

	print("running:", [exe_name])
	ret = subprocess.call([exe_name])
	if ret != 0:
		print("Failed to run executable %s, error %s"%(exe_name, ret))
		exit(ret)

	dt = numpy.dtype({"names":["t","p0","x","xp","y","yp","ct","dp"], "formats":["f8","f8","f8","f8","f8","f8","f8","f8"]})
	init_bunch = numpy.loadtxt("act_init_bunch.dat", skiprows=1, dtype=dt)

	dt = numpy.dtype({"names":["name","s","pos","x","xp","y","yp","type","id","loc", "turn"], "formats":["U32","f8","f8","f8","f8","f8","f8","i8","i8","f8","i8"]})
	losses = numpy.loadtxt("act_losses.dat", skiprows=1, dtype=dt)
	losses = numpy.sort(losses, order="id")

	if "plot" in sys.argv:
		from matplotlib import pyplot
		fig, ax = pyplot.subplots()
		ax.plot(init_bunch['x'], init_bunch['y'], "k.")
		cs = ["b","r","g","m","c","b"]
		for s in range(6):
			this_loss = losses[losses['s'] == s]
			ax.plot(this_loss['x'], this_loss['y'], cs[s]+"x")
		draw_rect(ax, "grey", 0.9, 0.8)
		draw_circ(ax, "grey", 0.75)
		draw_ellip(ax, "grey", 0.7, 0.6)
		draw_rectellip(ax, "grey", 0.4, 0.4, 0.4, 0.5)
		draw_octa(ax, "grey", 0.3, 0.3, pi/8, 3*pi/8)
		draw_octa(ax, "grey", 0.15, 0.2, pi/4, 3*pi/8)
		pyplot.show()

	loss_groups = [[len(list(g)),k] for k,g in itertools.groupby(losses["s"])]
	if "ref_run" in sys.argv:
		with open("aperture_config_test_ref.dat", "w") as outfile:
			print("# count, pos", file=outfile)
			for lg in loss_groups:
				print(lg[0], lg[1], file=outfile)

	dt = numpy.dtype({"names":["count","pos"], "formats":["i8","f8"]})
	loss_groups_ref = numpy.loadtxt(find_data_file("aperture_config_test_ref.dat"), skiprows=1, dtype=dt)

	good = True
	for n, lg in enumerate(loss_groups):
		if lg[0] != loss_groups_ref[n]["count"] or lg[1] != loss_groups_ref[n]["pos"]:
			print("Expected bin %s to have %s losses, got %s"%(n,lg,loss_groups_ref[n]))
			good = False
			break
	if not good:
		print("Failed")
		exit(1)
	else:
		print("Passed")
		os.remove("act_init_bunch.dat")
		os.remove("act_losses.dat")
		os.remove("act_aperture_conf.log")

def draw_rect(ax, m, hw, hh):
	ax.plot([hw,hw,-hw,-hw,hw], [hh,-hh,-hh,hh,hh], m)

def draw_circ(ax, m, r):
	angs = numpy.linspace(0,2*pi,200)
	xx, yy = r*numpy.sin(angs), r*numpy.cos(angs)
	ax.plot(xx,yy, m)

def draw_ellip(ax, m, ew, eh):
	angs = numpy.linspace(0,2*pi,200)
	xx, yy = ew*numpy.sin(angs), eh*numpy.cos(angs)
	ax.plot(xx,yy, m)

def draw_rectellip(ax, m, hw, hh, ew, eh):
	angs = numpy.linspace(0,2*pi,200)
	xx = (ew*numpy.sin(angs)).clip(-hw, hw)
	yy =  (eh*numpy.cos(angs)).clip(-hh, hh)
	ax.plot(xx, yy, m)

def draw_octa(ax, m, hw, hh, a1, a2):
	cy = hw * tan(a1)
	cx = hh * tan(pi/2 - a2)
	ax.plot([hw, hw, cx, -cx, -hw, -hw, -cx, cx, hw], [-cy, cy, hh, hh, cy, -cy, -hh, -hh, -cy], m)

if __name__ == "__main__":
	main()
