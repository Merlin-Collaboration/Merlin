#!/usr/bin/env python3

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2019 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING

import numpy
from matplotlib import pyplot
pyplot.rcParams['figure.figsize'] = 15, 10
from divnorm import DivergingNorm

import common

numpy.set_printoptions(linewidth=200)

def get_curves(tracker, element, settings=None):
	# x xp y yp ct dp
	pref = [0,0,0,0,0,0]
	particles = [pref.copy()]

	half_range = 1e-3
	nsteps = 21

	deltas = numpy.linspace(-half_range, half_range, nsteps)
	for c in range(6):
		for delta in deltas:
			p1 = pref.copy()
			p1[c] += delta
			particles.append(p1)
	
	p_out = tracker(element, particles, settings=settings)
	curves = numpy.zeros([6,6,nsteps,2], dtype="f8")

	p_in = numpy.array(particles)
	p_out = numpy.array(p_out)
	for i in range(6):
		for j in range(6):
			# initial j'th coord of particles with offsets in j'th
			curves[i,j,:,0] = p_in[1 + nsteps*j: 1 + nsteps*(j+1)][:,j]
			# final i'th coord of particles with offsets in j'th
			curves[i,j,:,1] = p_out[1 + nsteps*j: 1 + nsteps*(j+1)][:,i]
			
	return curves, particles, p_out


def plot_curves(trackers, name):
	fig, axs = pyplot.subplots(ncols=6, nrows=6, gridspec_kw={'hspace': 0.2, 'wspace': 0.2, "left":0.05, "right":0.95, "top":0.95, "bottom":0.05})

	colours = ["g", "r", "b", "c", "m"]
	for n, tracker in enumerate(trackers):
		axs[0,0].annotate(tracker["name"], [0.1+ 0.2*n, 0.98], xycoords="figure fraction", color=colours[n])

	for i in range(6):
		for j in range(6):
			for n, tracker in enumerate(trackers):
				cx = tracker["curves"][i,j,:,0]
				cy = tracker["curves"][i,j,:,1]
				x_off = cx.max() * 0.02 * n
				axs[i,j].plot(cx+x_off, cy, "-", c=colours[n])

			axs[i,j].annotate("[%s,%s]"%(i+1,j+1), [0.9,1], xycoords="axes fraction")
			axs[i,j].spines['bottom'].set_position(('data', 0))
			axs[i,j].spines['top'].set_visible(False)
			axs[i,j].spines['left'].set_position(('data', 0))
			axs[i,j].spines['right'].set_visible(False)

	out_name = "plots/"+name+".pdf"
	pyplot.savefig(out_name)
	print("Wrote", out_name)

def run_curves(run_settings, trackers, elements):
	for element in elements:
		el_name = element["type"] + "_" + "_".join("%s_%s"%(k,v) for k,v in sorted(element.items()) if k != "type")
		
		for t in trackers:
			full_settings = {}
			full_settings.update(run_settings)
			if "settings" in t:
				full_settings.update(t["settings"])
			curves, p_in, p_out = get_curves(t["f"], element, full_settings)

			t["curves"] = curves
			t["p_in"] = p_in
			t["p_out"] = p_out

		plot_curves(trackers, "multitrack_curves_"+el_name)
		continue
