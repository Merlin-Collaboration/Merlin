# Makefile for the MERLIN Examples
# -------------------------------------
# version 1.0 20.10.2005
# Nick Walker, DESY
#

SpinTracking : SpinTracking_Main.o 
	${CPP} ${LDFLAGS} -o $@ $^ -lmerlin -lm

%.o : %.cpp
	${CPP} ${CCFLAGS} -c -o $@ $<


