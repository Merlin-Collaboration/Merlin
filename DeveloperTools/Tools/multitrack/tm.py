#!/usr/bin/env python3

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2019 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING

import numpy
from matplotlib import pyplot
pyplot.rcParams['figure.figsize'] = 15, 7
from divnorm import DivergingNorm

import common

numpy.set_printoptions(linewidth=200)

def get_tm(tracker, element, settings=None):
	# x xp y yp ct dp
	pref = [0,0,0,0,0,0]
	particles = [pref.copy()]

	deltas = [1e-6] * 6
	for c in range(6):
		p1 = pref.copy()
		p2 = pref.copy()
		p1[c] -= deltas[c]
		p2[c] += deltas[c]
		
		particles.append(p1)
		particles.append(p2)
	
	p_out = tracker(element, particles, settings=settings)
	tm = numpy.eye(6, dtype="f8")

	for i in range(6):
		for j in range(6):
			if j == 4: continue# FIXME? madx_ptc ignores t input
			diff = p_out[1 + j*2 + 1][i] - p_out[1 + j*2][i]
			tm[i,j] = diff / 2 / deltas[j]
			
	return tm, particles, p_out

def plot_tm_diff(t1, t2, name):
	show_p0 = True
	
	if numpy.all(numpy.array([t1["p_out"][0], t2["p_out"][0]])==0):
		show_p0 =False
	
	if show_p0:
		fig, axs = pyplot.subplots(ncols=3, nrows=2,gridspec_kw={"height_ratios":[6,1]})
	else:
		fig, axs = pyplot.subplots(ncols=3)
	axs = axs.flatten()
	
	tm1, tm2 = t1["tm"], t2["tm"]
	tm_diff = numpy.abs(tm1 - tm2)
	divnorm = DivergingNorm(vmin=-1, vcenter=0, vmax=1)
	divnorm2 = DivergingNorm(vmin=-1e-3, vcenter=0, vmax=1e-3)
	axs[0].imshow(tm1, interpolation="nearest", cmap="seismic", norm=divnorm)
	axs[1].imshow(tm2, interpolation="nearest", cmap="seismic", norm=divnorm)
	axs[2].imshow(tm_diff, interpolation="nearest", cmap="seismic", norm=divnorm2)

	axs[0].set_title(t1["name"])
	axs[1].set_title(t2["name"])

	for i in range(6):
		for j in range(6):
			axs[0].annotate("[%s,%s]"%(i+1,j+1), [j-0.3,i-0.2], backgroundcolor="white")
			axs[1].annotate("[%s,%s]"%(i+1,j+1), [j-0.3,i-0.2], backgroundcolor="white")
			axs[2].annotate("[%s,%s]"%(i+1,j+1), [j-0.3,i-0.2], backgroundcolor="white")

			axs[0].annotate("%.3g"%(tm1[i,j]), [j-0.3,i+0.1], backgroundcolor="white")
			axs[1].annotate("%.3g"%(tm2[i,j]), [j-0.3,i+0.1], backgroundcolor="white")
			axs[2].annotate("%.3g"%(tm_diff[i,j]), [j-0.3,i+0.1], backgroundcolor="white")

	if show_p0:
		p0_1, p0_2 = t1["p_out"][0], t2["p_out"][0]
		p0_1 = numpy.array(p0_1, ndmin=2)
		p0_2 = numpy.array(p0_2, ndmin=2)
		p0_diff = numpy.abs(p0_1 - p0_2)
		axs[3].imshow(p0_1, interpolation="nearest", cmap="seismic", norm=divnorm)
		axs[4].imshow(p0_2, interpolation="nearest", cmap="seismic", norm=divnorm)
		axs[5].imshow(p0_diff, interpolation="nearest", cmap="seismic", norm=divnorm2)
		for i in range(6):
			axs[3].annotate("%.3g"%(p0_1[0,i]),    [i-0.3,0.1], backgroundcolor="white")
			axs[4].annotate("%.3g"%(p0_2[0,i]),    [i-0.3,0.1], backgroundcolor="white")
			axs[5].annotate("%.3g"%(p0_diff[0,i]), [i-0.3,0.1], backgroundcolor="white")
		
		for ax in axs[3:6]:
			ax.yaxis.set_ticks([])
			ax.xaxis.set_ticks([0,1,2,3,4,5])
			ax.xaxis.set_ticklabels(["x", "xp", "y", "yp", "ct", "dp"])

	pyplot.tight_layout()
	out_name = "plots/"+name+".pdf"
	pyplot.savefig(out_name)
	print("Wrote", out_name)
	#pyplot.show()

def plot_tms(tms, name):
	show_p0 = True
	ncols = len(tms)
	
	if show_p0:
		fig, axs = pyplot.subplots(ncols=ncols, nrows=2, gridspec_kw={"height_ratios":[6,1]})
	else:
		fig, axs = pyplot.subplots(ncols=ncols)
	axs = axs.flatten()
	
	divnorm = DivergingNorm(vmin=-1, vcenter=0, vmax=1)
	divnorm2 = DivergingNorm(vmin=-1e-3, vcenter=0, vmax=1e-3)
	
	for n, tm in enumerate(tms):
		axs[n].imshow(tm["tm"], interpolation="nearest", cmap="seismic", norm=divnorm)
		axs[n].set_title(tm["name"])

		for i in range(6):
			for j in range(6):
				axs[n].annotate("[%s,%s]"%(i+1,j+1), [j-0.3,i-0.2], backgroundcolor="white")
				axs[n].annotate("%.3g"%(tm["tm"][i,j]), [j-0.3,i+0.1], backgroundcolor="white")

		if show_p0:
			axs[ncols+n].imshow(numpy.array(tm["p_out"][0], ndmin=2), interpolation="nearest", cmap="seismic", norm=divnorm)
			for i in range(6):
				axs[ncols+n].annotate("%.3g"%(tm["p_out"][0][i]),  [i-0.3,0.1], backgroundcolor="white")

			axs[ncols+n].yaxis.set_ticks([])
			axs[ncols+n].xaxis.set_ticks([0,1,2,3,4,5])
			axs[ncols+n].xaxis.set_ticklabels(["x", "xp", "y", "yp", "ct", "dp"])

	pyplot.tight_layout()
	out_name = "plots/"+name+".pdf"
	pyplot.savefig(out_name)
	print("Wrote", out_name)

def run_transfermatrix(run_settings, trackers, elements):
	for element in elements:
		el_name = element["type"] + "_" + "_".join("%s_%s"%(k,v) for k,v in sorted(element.items()) if k != "type")
		
		for t in trackers:
			full_settings = {}
			full_settings.update(run_settings)
			if "settings" in t:
				full_settings.update(t["settings"])
			tm, p_in, p_out = get_tm(t["f"], element, full_settings)

			t["tm"] = tm
			t["p_in"] = p_in
			t["p_out"] = p_out

		for t in trackers:
			print(t["name"])
			print(t["tm"])

		if len(trackers) == 2:
			plot_tm_diff(trackers[0], trackers[1], "multitrack_tm_diff_"+el_name)
		else:
			plot_tms(trackers, "multitrack_tm_"+el_name)
		continue
