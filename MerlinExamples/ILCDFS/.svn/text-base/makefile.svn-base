# Makefile for the MERLIN Examples
# -------------------------------------
# version ILCDFS, 14.1.08
# Dirk Kruecker, DESY
#
VPATH = src


OFILES = ILCDFS_main.o Accelerator.o AcceleratorWithErrors.o \
         ConstantGradientAdjustment.o ConstructSrot.o DFSApp.o \
         DFSCorrection.o DFSOutput.o EnergyAdjustmentPolicy.o \
         ILCDFS_IO.o ModelConstruction.o OneToOneCorrection.o \
         OptimumDispersionMatch.o ParticleTrackingModel.o \
         ResponseMatrixGenerator.o SMPTrackingModel.o \
         TeslaWakePotential.o XTFFInterface_1.o


ILCDFS : ${OFILES}
	${CPP} ${LDFLAGS} -o $@ $^ -lmerlin -lm

%.o : %.cpp
	${CPP} ${CCFLAGS} -c -o $@ $<

