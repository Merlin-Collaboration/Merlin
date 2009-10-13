# Makefile for the MERLIN Examples
# -------------------------------------
# version 1.0 20.12.2004
# Nick Walker, DESY
#

OFILES = HaloCollimation_main.o HaloTracker.o QuadIntegrator.o ConstructModel.o

HaloCollimation : ${OFILES}
	${CPP} ${LDFLAGS} -o $@ $^ -lmerlin -lm

%.o : %.cpp
	${CPP} ${CCFLAGS} -c -o $@ $<


