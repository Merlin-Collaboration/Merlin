# Makefile for the MERLIN Examples
# -------------------------------------
# version 1.0 20.12.2004
# Nick Walker, DESY
#

OFILES = model_construction.o  \
	 QuadReferenceOutput.o \
	 TeslaWakePotential.o  \
	 TrackingOutput.o      \
         ILCML_main.o

ILCML : ${OFILES}
	${CPP} ${LDFLAGS} -o $@ $^ -lmerlin -lm

%.o : %.cpp
	${CPP} ${CCFLAGS} -c -o $@ $<


