#!/usr/bin/env python3

# Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
# Copyright (c) 2001-2019 The Merlin++ developers
# This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING

import tempfile
import os
import errno

this_dir = os.path.dirname(__file__)

tmp_dir = tempfile.TemporaryDirectory(prefix="multitrack_")
tmp_dir_path = tmp_dir.name
print("Temp dir:", tmp_dir_path)

def close_temp():
	tmp_dir.cleanup()

import atexit
atexit.register(close_temp)

settings = {
	"energy": 7000,
	"particle": "proton",
}

quiet = False

def mkdir(path):
	try:
		os.makedirs(path)
	except OSError as e:
		if e.errno == errno.EEXIST:
			pass
		else:
			raise
