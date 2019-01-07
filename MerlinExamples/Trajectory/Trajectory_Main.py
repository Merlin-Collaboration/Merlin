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

sim_filename = "Trajectory.dat"

for run_n in range(1):
    # multiple runs useful for visualising if differences are statistical or systematic
    
    # Run the test script
    try:
        os.remove(sim_filename)
    except OSError:
        pass
    print("Running Trajectory")
    exe_name = os.path.join(os.path.dirname(__file__), "Trajectory")
    print("running:", [exe_name])
    ret = subprocess.call([exe_name])
        
    if ret != 0:
        print("Failed to run executable %s, error %s"%(exe_name, ret))
        exit(ret)
    
    try:
        tracking_data = numpy.loadtxt(sim_filename)
    except IOError:
        print("Failed to read", sim_filename)
        exit(1)
    print("Read simulated data")
    
    if "plot" in sys.argv:
        from matplotlib import pyplot
        fig, axes = plt.subplots(nrows=2, sharex=False)
        fig.set_size_inches(8,11)
        axes[0].set_title('X')
        axes[0].plot(data[:,0])
        axes[1].set_title('Y')
        axes[1].plot(data[:,2]) 
        
if "plot" in sys.argv:
    pyplot.show()

if all_good:
    print("Pass")
    os.remove(sim_filename)
else:
    print("Fail")
    exit(1)
        