#!/usr/bin/env python3
import sys

import numpy
from matplotlib import pyplot

cols = "T P0 X XP Y YP CT DP".split()

def read_bunch(fname):
	bunch = numpy.genfromtxt(fname, dtype=[(s, "f8") for s in cols])
	return bunch



bunch1 = read_bunch("example_transferline_bunch_inject.out")
bunch2 = read_bunch("example_transferline_bunch_end.out")


pyplot.plot(bunch1['X'], bunch1['XP'], "b.", label="inject")
pyplot.plot(bunch2['X'], bunch2['XP'], "r.", label="end")
pyplot.legend()
pyplot.show()
