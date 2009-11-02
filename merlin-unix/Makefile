# the compiler to use
export CPP= g++4
#export CPP=distcc
# the root directory for the Merlin source code
export MERLIN_ROOT_DIR=${PWD}/../Merlin

# the root build directory
# [see README for more details]
export MERLIN_BUILD_DIR=${PWD}

# the install directory where the libs will be places
export MERLIN_INSTALL_DIR=${HOME}/lib

export INCPATH= -I${MERLIN_ROOT_DIR}
export CCFLAGS= -pipe -fPIC

export OBJDIR=${MERLIN_BUILD_DIR}/obj/release
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

default: release

.PHONY: subdirs $(SUBDIRS)
subdirs: $(SUBDIRS)
$(SUBDIRS): obj
	$(MAKE) -C ${MERLIN_ROOT_DIR}/$@

obj: config
	if [ ! -d ${MERLIN_BUILD_DIR}/obj ] ; \
	then \
		mkdir ${MERLIN_BUILD_DIR}/obj ; \
	fi
	if [ ! -d ${MERLIN_BUILD_DIR}/obj/release ] ; \
	then \
		mkdir ${MERLIN_BUILD_DIR}/obj/release ; \
	fi

.PHONY: release
release: subdirs
	gcc -shared -Wl,--export-dynamic,-soname=libmerlin.so -o ${MERLIN_INSTALL_DIR}/libmerlin.so ${MERLIN_BUILD_DIR}/obj/release/*.o

clean:
	if [ -d ${MERLIN_BUILD_DIR}/obj ]; then rm -r ${MERLIN_BUILD_DIR}/obj; fi

	for dirn in ${SUBDIRS}; do \
		cd ${MERLIN_ROOT_DIR}/$${dirn}; \
		rm -f *~ catagory Makefile; \
	done

config:
	for dirn in ${SUBDIRS}; do \
		cd ${MERLIN_ROOT_DIR}/$${dirn}; \
		echo "building makefile for " $${dirn}; \
		${MERLIN_BUILD_DIR}/build_makefile; \
	done
