#!/usr/bin/env python3

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2019 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING

import subprocess
import os

import common
from madx import make_element
from merlin import make_particle

merlin_command =  os.path.join(common.this_dir, "merlin_track")
madinterface_tfs_path = "madinterface.tfs"

madx_temp = """
option, -echo,  warn, -info, verify;
M1: MARKER;
M2: MARKER;
TEST_ELEMENT: {element};

TEST_CELL: LINE=(M1, TEST_ELEMENT, M2);

BEAM, PARTICLE={particle}, ENERGY={energy}, SEQUENCE=TEST_CELL;
USE, PERIOD = TEST_CELL;

//SAVE, FILE=save, SEQUENCE=TEST_CELL;
//USE, PERIOD = TEST_CELL; //  need to recall USE after SAVE


//DUMPSEQU , SEQUENCE=TEST_CELL, LEVEL=3;
//stop;

//SELECT, flag=survey, //column=NAME,KEYWORD,S,L,KS,KSL,K0L,K1L,K2L,K3L,K4L,K1S,K2S,K3S,K4S,HKICK,VKICK,BETX,BETY,ALFX,ALFY,MUX,MUY,DX,DY,DPX,DPY,R11,R12,R22,R21,X,PX,Y,PY,T,DELTAP,VOLT,LAG,HARMON,FREQ,E1,E2,APERTYPE,APER_1,APER_2,APER_3,APER_4,TILT,ANGLE,assembly_ID,mech_sep, samsam;
SURVEY, FILE={madinterface_tfs_path};

//stop;

SELECT, flag=twiss,clear;
SELECT, flag=twiss, column=NAME,KEYWORD,S,L,KS,KSL,K0L,K1L,K2L,K3L,K4L,K1S,K2S,K3S,K4S,HKICK,VKICK,BETX,BETY,ALFX,ALFY,MUX,MUY,DX,DY,DPX,DPY,R11,R12,R22,R21,X,PX,Y,PY,T,DELTAP,VOLT,LAG,HARMON,FREQ,E1,E2,APERTYPE,APER_1,APER_2,APER_3,APER_4,TILT,ANGLE,assembly_ID,mech_sep;
TWISS, sequence=TEST_CELL, centre, file={madinterface_tfs_path}, BETX=1, BETY=1;


stop;

"""

def make_mad_tfs(element, settings=None):
	new_element = {}
	new_element.update(null_element)
	new_element.update(element)
	element = new_element
	run_settings = common.settings.copy()
	if settings is not None:
		run_settings.update(settings)

	mad_in = madx_temp.format(element=make_element(element),
	                          madinterface_tfs_path=madinterface_tfs_path,
	                          **run_settings)
	if common.quiet:
		stdout = subprocess.PIPE
	else:
		stdout = None

	proc = subprocess.Popen(settings["madx_command"], stdin=subprocess.PIPE, stdout=stdout, universal_newlines=True, cwd=common.tmp_dir_path)
	proc.communicate(input=mad_in)
	ret = proc.wait()
	if ret != 0:
		raise Exception("Madx failed")
	
	for line in open(os.path.join(common.tmp_dir_path, madinterface_tfs_path)):
		if "TEST_ELEMENT" in line:
			print(line)



null_element={"Len":0, "k1":0, "k2":0}
def track_particles(element, particles, settings=None):
	ps = []
	new_element = {}
	new_element.update(null_element)
	new_element.update(element)
	element = new_element

	make_mad_tfs(element, settings)

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
		merlin_infile.write("madinterface %s"%madinterface_tfs_path)
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
	os.remove(os.path.join(common.tmp_dir_path, madinterface_tfs_path))
	return outbunch
