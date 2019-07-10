#!/usr/bin/env python3

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2019 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING

import subprocess
import os

import common

madx_temp = """
option, -echo,  warn, -info, verify;
M1: MARKER;
M2: MARKER;
TEST_ELEMENT: {element};

TEST_CELL: LINE=(M1, TEST_ELEMENT, M2);

BEAM, PARTICLE={particle}, ENERGY={energy}, SEQUENCE=TEST_CELL;
USE, PERIOD = TEST_CELL;

{mad_track}

stop;
"""

madx_track_ptc = """
ptc_create_universe;
ptc_create_layout,model=1,method=2,nst=1,exact=false;

{beam}

ptc_track,icase=6,closed_orbit,dump, turns=1,ffile=1, ELEMENT_BY_ELEMENT=true, onetable, CLOSED_ORBIT=false;

ptc_track_end;
ptc_end;
system, "grep -v '^#segment' trackone > ptc_track.tfs"

"""

madx_track_thin = """
SELECT, CLASS=HKICKER, FLAG=makethin, SLICE=2;
MAKETHIN, SEQUENCE=TEST_CELL;
TRACK, ONEPASS=true, ONETABLE=true, FILE=trackthin;

{beam}

RUN, TURNS=1;

ENDTRACK;
system, "grep -v '^#segment' trackthinone > thin_track.tfs"

"""


def make_element(element):
	print("element", element)
	if element["type"] == "drift":
		out = "DRIFT, L=%s" % (element["len"])
	elif element["type"] == "quad":
		out = "QUADRUPOLE, L={len}, K1={k1}".format(**element)
	elif element["type"] == "sex":
		out = "SEXTUPOLE, L={len}, K2={k2}".format(**element)
	elif element["type"] == "oct":
		out = "OCTUPOLE, L={len}, K3={k3}".format(**element)
	elif element["type"] == "sbend":
		out = "SBEND, L={len}, ANGLE={angle}, K1={k1}, K2={k2}".format(**element)
	elif element["type"] == "vkick":
		out = "VKICKER, L={len}, KICK={vkick}".format(**element)
	elif element["type"] == "hkick":
		out = "HKICKER, L={len}, KICK={hkick}".format(**element)
	elif element["type"] == "multipole":
		knl_s = ", ".join([str(x) for x in element["knl"]])
		ksl_s = ", ".join([str(x) for x in element["ksl"]])
		out = "MULTIPOLE, KNL={{{knl_s}}}, KSL={{{ksl_s}}}".format(**element, knl_s=knl_s, ksl_s=ksl_s)
	else:
		raise Exception("Unknown element: %s"%element["type"])
	return out

def make_particle_thin(particle):
	out = "START, x={}, px={}, y={}, py={}, t={}, pt={};".format(*particle)
	return out

def make_particle_ptc(particle):
	out = "ptc_start, x={}, px={}, y={}, py={}, t={}, pt={};".format(*particle)
	return out

null_element={"Len":0, "k1":0, "k2":0}
def track_particles(element, particles, settings=None):
	ps = []
	new_element = {}
	new_element.update(null_element)
	new_element.update(element)
	element = new_element

	run_settings = common.settings.copy()
	if settings is not None:
		run_settings.update(settings)

	if settings["mad_mode"] == "ptc":
		for p in particles:
			ps.append(make_particle_ptc(p))
		beam = "\n".join(ps)
		mad_track = madx_track_ptc.format(beam=beam, **run_settings)
		mad_in = madx_temp.format(element=make_element(element), mad_track=mad_track,
	                          **run_settings)
		madx_outfile_path = os.path.join(common.tmp_dir_path,"ptc_track.tfs")
		madx_tmp_path = os.path.join(common.tmp_dir_path,"trackone")
	elif settings["mad_mode"] == "thin":
		for p in particles:
			ps.append(make_particle_thin(p))
		beam = "\n".join(ps)
		mad_track = madx_track_thin.format(beam=beam, **run_settings)
		mad_in = madx_temp.format(element=make_element(element), mad_track=mad_track,
	                          **run_settings)
		madx_outfile_path = os.path.join(common.tmp_dir_path,"thin_track.tfs")
		madx_tmp_path = os.path.join(common.tmp_dir_path,"trackthinone")
	else:
		raise ValueError("Unknown mad_mode for madx:", settings["mad_mode"])

	if common.quiet:
		stdout = subprocess.PIPE
	else:
		stdout = None

	proc = subprocess.Popen(settings["madx_command"], stdin=subprocess.PIPE, stdout=stdout, universal_newlines=True, cwd=common.tmp_dir_path)
	proc.communicate(input=mad_in)
	ret = proc.wait()
	if ret != 0:
		raise Exception("Madx failed")

	outbunch = []
	with open(madx_outfile_path) as track_file:
		for line in track_file:
			if line[0] in "@*$": continue
			number, turn, x, px, y, py, t, pt, s, e = [float(_x) for _x in line.split()]
			
			if int(turn)>0:
				outbunch.append([x, px, y, py, t, pt])
	if len(particles) != len(outbunch):
		raise Exception("Particles lost: start: %s end: %s"%(len(particles), len(outbunch)))

	os.remove(madx_outfile_path)
	os.remove(madx_tmp_path)
	return outbunch
