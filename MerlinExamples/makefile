# Makefile for the MERLIN Examples
# -------------------------------------
# version 2.0
# DK 14.1.2009
# version 1.1 20.10.2005
# Nick Walker, DESY
#

#----- Edit these lines ---------------------------------------
# the compiler to use
CPP= g++
# the root directory for the Merlin source code
MERLIN_ROOT_DIR=${PWD}/../Merlin
MERLIN_LIB_PATH=${HOME}/lib
#---------------------------------------------------------------

INCPATH= -I${MERLIN_ROOT_DIR}
CCFLAGS= -pg -march=native -mfpmath=sse -O3 -fopenmp ${INCPATH}
LDFLAGS= -L${MERLIN_LIB_PATH} -lgsl -lgslcblas -fopenmp

EGDIRS= ClosedOrbit LatticeFunctions StrayFieldClosedOrbit \
        Tracking Trajectory HaloCollimation SpinTracking ILCML ILCDFS \
        ROOT_trackingOutput Wakefields AcceleratorErrors

merlin_eg:
	for dirn in ${EGDIRS}; do \
		cd $${dirn}; \
		echo "building catagory " $${dirn}; \
		${MAKE}  "CPP=${CPP}" "CCFLAGS=${CCFLAGS}" \
	                  "LDFLAGS=${LDFLAGS}" ; \
                cd .. ; \
	done

clean:
	for dirn in ${EGDIRS}; do \
		cd $${dirn}; \
		rm -f *.o *~ *.log *.out; \
		cd .. ; \
	done
