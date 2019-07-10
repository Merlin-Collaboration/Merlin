#!/usr/bin/env python3

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2019 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING

import subprocess
import os

import common

merlin_command =  os.path.join(common.this_dir, "merlin_track")

def make_particle(particle):
	out = "particle {} {} {} {} {} {}".format(*particle)
	return out

def make_element(element):
	if element["type"] == "drift":
		out = "drift %s" % (element["len"])
	elif element["type"] == "quad":
		out = "quad {len} {k1}".format(**element)
	elif element["type"] == "sbend":
		out = "sbend {len} {angle} {k1} {k2}".format(**element)
	elif element["type"] == "vkick":
		out = "vkick {len} {vkick}".format(**element)
	elif element["type"] == "hkick":
		out = "hkick {len} {hkick}".format(**element)
	else:
		raise Exception("Unknown element: %s"%element["type"])
	return out

null_element={"Len":0, "k1":0, "k2":0}
def track_particles(element, particles, settings=None):
	ps = []
	new_element = {}
	new_element.update(null_element)
	new_element.update(element)
	element = new_element

	for p in particles:
		ps.append(make_particle(p))
	beam = "\n".join(ps)

	run_settings = common.settings.copy()
	if settings is not None:
		run_settings.update(settings)

	merlin_infile_path = os.path.join(common.tmp_dir_path, "merlin_in.dat")
	with open(merlin_infile_path, "w") as merlin_infile:
		for k,v in run_settings.items():
			merlin_infile.write("set %s %s\n"%(k,v))
		merlin_infile.write(beam)
		merlin_infile.write("\n")
		merlin_infile.write(make_element(element))
		merlin_infile.write("\n")

	command = [merlin_command, merlin_infile_path]
	print("Running", command)
	proc = subprocess.Popen(command,
	                        universal_newlines=True, cwd=common.tmp_dir_path)

	ret = proc.wait()
	if ret != 0:
		raise Exception("Merlin failed")

	outbunch = []
	merlin_outfile_path = os.path.join(common.tmp_dir_path,"merlin_out.dat")
	with open(merlin_outfile_path) as track_file:

		for line in track_file:
			if line[0] in "@*$#": continue
			if line.strip() == "": continue
			x, px, y, py, t, pt = [float(_x) for _x in line.split()]

			outbunch.append([x, px, y, py, t, pt])

	if len(particles) != len(outbunch):
		raise Exception("Particles lost: start: %s end: %s"%(len(particles), len(outbunch)))

	os.remove(merlin_infile_path)
	os.remove(merlin_outfile_path)
	return outbunch
