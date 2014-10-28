# the compiler to use
#export CPP= g++451
export CPP= g++
#export CPP=distcc
#export LINK=gcc451
export LINK=gcc
#export LINK=gcc45svn

# the root directory for the Merlin source code
export MERLIN_ROOT_DIR=${PWD}/../Merlin

# the root build directory
# [see README for more details]
export MERLIN_BUILD_DIR=${PWD}

# the install directory where the libs will be placed
export MERLIN_INSTALL_DIR=${HOME}/lib

export INCPATH= -I${MERLIN_ROOT_DIR} -I${HOME}/include

#export CCFLAGS=-pipe -fPIC -O2 -march=nocona -msse -msse2 -mfpmath=sse -pedantic -Wall -lgsl -L${HOME}/lib

export CCFLAGS=-flto -pipe -fPIC -O3 -march=native -ffast-math -std=c++11 -msse -msse2 -mfpmath=sse,387 -pedantic -Wall #-lgsl -L${HOME}/lib
#export CCFLAGS=-pipe -fPIC -ggdb -g -pg -pedantic -Wall -lgsl -L${HOME}/lib

#export CCFLAGS=-pipe -fPIC -pedantic -Wall #-lgsl -L${HOME}/lib 
export OBJDIR=${MERLIN_BUILD_DIR}/obj/

SUBDIRS= AcceleratorModel \
      AcceleratorModel/ActiveMonitors \
      AcceleratorModel/Apertures \
      AcceleratorModel/Construction \
      AcceleratorModel/ControlElements \
      AcceleratorModel/Frames \
      AcceleratorModel/Implementation \
      AcceleratorModel/Implementation/BPMChannels \
      AcceleratorModel/Miscellaneous \
      AcceleratorModel/StdComponent \
      AcceleratorModel/StdField \
      AcceleratorModel/StdGeometry \
      AcceleratorModel/Supports \
      AcceleratorModel/TrackingInterface \
      BasicTransport \
      BeamDynamics \
      BeamDynamics/ParticleTracking \
      BeamDynamics/ParticleTracking/Integrators \
      BeamDynamics/SMPTracking \
      BeamDynamics/CommonUtilities \
      BeamModel \
      Channels \
      Collimators \
      TLAS \
      Corrections \
      EuclideanGeometry \
      GroundMotionModels \
      IO \
      MADInterface \
      NumericalUtils \
      Random \
      utility \
      RingDynamics 

ifeq ($(build),debug)
	CCFLAGS:=${CCFLAGS} -ggdb3 -ggdb -g3 -W -Wextra -fopenmp -DENABLE_OPENMP
	OBJDIR:= $(OBJDIR)debug
else ifeq ($(build),release)
	CCFLAGS:=${CCFLAGS} -DNDEBUG -fopenmp -DENABLE_OPENMP -O3 -march=native -mfpmath=sse,387
	OBJDIR:= $(OBJDIR)release
else ifeq ($(build),mpi)
	##set to the compiler one wishes openMPI to use
	export OMPI_CXX=g++
	export OMPI_CC=gcc
	CPP=mpiCC
	LINK=mpiCC
	##CPP=mpicc
	##LINK=mpicc
	CCFLAGS:= ${CCFLAGS} -DENABLE_MPI -DNDEBUG -O3 -march=native -mfpmath=sse,387
	#CCFLAGS:= ${CCFLAGS} -DENABLE_MPI
	OBJDIR:= $(OBJDIR)mpi
else
	##Default to release if no type specified
	CCFLAGS:=${CCFLAGS} -mfpmath=sse -DNDEBUG -fopenmp -DENABLE_OPENMP
	OBJDIR:= $(OBJDIR)release
endif
default: build

.PHONY: subdirs $(SUBDIRS)
subdirs: $(SUBDIRS)
$(SUBDIRS): obj
	$(MAKE) -C ${MERLIN_ROOT_DIR}/$@

obj:
	if [ ! -d ${MERLIN_BUILD_DIR}/obj ] ; \
	then mkdir ${MERLIN_BUILD_DIR}/obj ; \
	fi
	if [ ! -d ${MERLIN_BUILD_DIR}/obj/release ] ; \
	then mkdir ${MERLIN_BUILD_DIR}/obj/release ; \
	fi
	if [ ! -d ${MERLIN_BUILD_DIR}/obj/debug ] ; \
	then mkdir ${MERLIN_BUILD_DIR}/obj/debug; \
	fi
	if [ ! -d ${MERLIN_BUILD_DIR}/obj/mpi ] ; \
	then mkdir ${MERLIN_BUILD_DIR}/obj/mpi; \
	fi

build: subdirs
	$(LINK) -flto -fPIC -shared -Wl,--export-dynamic,-soname=libmerlin.so -o ${MERLIN_INSTALL_DIR}/libmerlin.so ${OBJDIR}/*.o
# -lrt
#	uncomment the following line for static linking
#	ar rcsv ${MERLIN_INSTALL_DIR}/libmerlin.a ${OBJDIR}/*.o
release:
	$(MAKE) $(MAKEFILE) build=release

debug:
	$(MAKE) $(MAKEFILE) build=debug

mpi:
	$(MAKE) $(MAKEFILE) build=mpi

.PHONY: clean
clean:
	if [ -d ${MERLIN_BUILD_DIR}/obj ]; then rm -r ${MERLIN_BUILD_DIR}/obj; fi

	for dirn in ${SUBDIRS}; do \
		cd ${MERLIN_ROOT_DIR}/$${dirn}; \
		rm -f *~ catagory Makefile; \
	done

.PHONY: config
config:
	for dirn in ${SUBDIRS}; do \
		cd ${MERLIN_ROOT_DIR}/$${dirn}; \
		echo "building makefile for " $${dirn}; \
		sh ${MERLIN_BUILD_DIR}/build_makefile; \
	done
