# Makefile for the MERLIN Examples
# -------------------------------------
# version 1.0 20.12.2004
# Nick Walker, DESY
#

StrayFieldClosedOrbit : StrayFieldClosedOrbit_Main.o ConstantStrayFieldProcess.o
	${CPP} ${LDFLAGS} -o $@ $^ -lmerlin -lm

%.o : %.cpp
	${CPP} ${CCFLAGS} -c -o $@ $<


