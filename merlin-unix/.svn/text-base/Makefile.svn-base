# Makefile for the MERLIN class library
# -------------------------------------
# Modified 25.11.08
# Dirk Kruecker
#
# version 1.2 22.11.01
# Nick Walker, DESY
#
# This is a simple makefile for configuring and compiling the Merlin
# source code. You should customise the indicated lines to suite your
# requirements
#
# This file has been tested using GNU's gmake. 
#
# See README file for more details
#
#----- Edit these lines -----------

# the compiler to use
CPP= g++

# the root directory for the Merlin source code
MERLIN_ROOT_DIR=${PWD}/../Merlin

# the root build directory
# [see README for more details]
MERLIN_BUILD_DIR=${PWD}

# the install directory where the libs will be places
MERLIN_INSTALL_DIR=${MERLIN_BUILD_DIR}/lib

#----------------------------------

INCPATH= -I${MERLIN_ROOT_DIR} 
#CCFLAGS= -Wall -malign-double -ffast-math -march=pentiumpro
#CCFLAGS= -w -malign-double -ffast-math -march=pentiumpro
#CCFLAGS= -w -fpermissive
CCFLAGS= -w 

CATS= AcceleratorModel \
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
      BeamDynamics/ParticleTracking \
      BeamDynamics/SMPTracking \
      BeamDynamics/CommonUtilities \
      BeamModel \
      Channels \
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
#      LinearOptics \

CATS1 = TLAS

dummy:

install: config depend release

install_debug: config depend debug

clean:
	if [ -d ${MERLIN_BUILD_DIR}/obj ]; then rm -r ${MERLIN_BUILD_DIR}/obj; fi

config:
	for dirn in ${CATS}; do \
		cd ${MERLIN_ROOT_DIR}/$${dirn}; \
		echo "building makefile for " $${dirn}; \
		${MERLIN_BUILD_DIR}/build_makefile; \
	done

depend:
	for dirn in ${CATS}; do \
		cd ${MERLIN_ROOT_DIR}/$${dirn}; \
		echo "updating dependencies for " $${dirn}; \
		${MERLIN_BUILD_DIR}/make_depend ${MERLIN_ROOT_DIR}; \
	done


debug:
	if [ ! -d ${MERLIN_BUILD_DIR}/obj ] ; \
	then \
		mkdir ${MERLIN_BUILD_DIR}/obj ; \
	fi
	if [ ! -d ${MERLIN_BUILD_DIR}/obj/debug ] ; \
	then \
		mkdir ${MERLIN_BUILD_DIR}/obj/debug ; \
	fi
	OBJDIR=${MERLIN_BUILD_DIR}/obj/debug
	for dirn in ${CATS}; do \
		cd ${MERLIN_ROOT_DIR}/$${dirn}; \
		echo "building catagory " $${dirn}; \
		${MAKE} catagory "CPP=${CPP}" "INCPATH=${INCPATH}" \
			"CCFLAGS=${CCFLAGS} -g2" \
			"OBJDIR=${MERLIN_BUILD_DIR}/obj/debug" ; \
	done
	if [ ! -d ${MERLIN_BUILD_DIR}/lib ] ; then mkdir ${MERLIN_BUILD_DIR}/lib ; fi
	if [ -f ${MERLIN_BUILD_DIR}/lib/libmerlin_debug.a ] ; then rm ${MERLIN_BUILD_DIR}/lib/libmerlin_debug.a; fi
	ar rv ${MERLIN_BUILD_DIR}/lib/libmerlin_debug.a ${MERLIN_BUILD_DIR}/obj/debug/*.o
	# move to install directory
	if [ ${MERLIN_INSTALL_DIR}/libmerlin_debug.a != ${MERLIN_BUILD_DIR}/lib/libmerlin_debug.a ] ; then \
		mv ${MERLIN_BUILD_DIR}/lib/libmerlin_debug.a ${MERLIN_INSTALL_DIR}/libmerlin_debug.a ; \
	fi
release:
	${CPP} -v;
	if [ ! -d ${MERLIN_BUILD_DIR}/obj ] ; \
	then \
		mkdir ${MERLIN_BUILD_DIR}/obj ; \
	fi
	if [ ! -d ${MERLIN_BUILD_DIR}/obj/release ] ; \
	then \
		mkdir ${MERLIN_BUILD_DIR}/obj/release ; \
	fi
	OBJDIR=${MERLIN_BUILD_DIR}/obj/release
	for dirn in ${CATS}; do \
		cd ${MERLIN_ROOT_DIR}/$${dirn}; \
		echo "building catagory " $${dirn}; \
		${MAKE} catagory "CPP=${CPP}" "INCPATH=${INCPATH}" \
			"CCFLAGS=${CCFLAGS} -O3 -DNDEBUG" \
			"OBJDIR=${MERLIN_BUILD_DIR}/obj/release" ; \
	done
	if [ ! -d ${MERLIN_BUILD_DIR}/lib ] ; then mkdir ${MERLIN_BUILD_DIR}/lib ; fi
	if [ -f ${MERLIN_BUILD_DIR}/lib/libmerlin.a ] ; then rm ${MERLIN_BUILD_DIR}/lib/libmerlin.a; fi
	ar rv ${MERLIN_BUILD_DIR}/lib/libmerlin.a ${MERLIN_BUILD_DIR}/obj/release/*.o
	# move to install directory
	if [ ${MERLIN_INSTALL_DIR}/libmerlin.a != ${MERLIN_BUILD_DIR}/lib/libmerlin.a ] ; then \
		mv ${MERLIN_BUILD_DIR}/lib/libmerlin.a ${MERLIN_INSTALL_DIR}/libmerlin.a ; \
	fi
	
profile:
	if [ ! -d ${MERLIN_BUILD_DIR}/obj ] ; \
	then \
		mkdir ${MERLIN_BUILD_DIR}/obj ; \
	fi
	if [ ! -d ${MERLIN_BUILD_DIR}/obj/profile ] ; \
	then \
		mkdir ${MERLIN_BUILD_DIR}/obj/profile ; \
	fi
	OBJDIR=${MERLIN_BUILD_DIR}/obj/profile
	for dirn in ${CATS}; do \
		cd ${MERLIN_ROOT_DIR}/$${dirn}; \
		echo "building catagory " $${dirn}; \
		${MAKE} catagory "CPP=${CPP}" "INCPATH=${INCPATH}" \
			"CCFLAGS=${CCFLAGS} -O3 -pg" \
			"OBJDIR=${MERLIN_BUILD_DIR}/obj/profile" ; \
	done
	if [ ! -d ${MERLIN_BUILD_DIR}/lib ] ; then mkdir ${MERLIN_BUILD_DIR}/lib ; fi
	if [ -f ${MERLIN_BUILD_DIR}/lib/libmerlin_pg.a ] ; then rm ${MERLIN_BUILD_DIR}/lib/libmerlin_pg.a; fi
	ar rv ${MERLIN_BUILD_DIR}/lib/libmerlin_pg.a ${MERLIN_BUILD_DIR}/obj/profile/*.o
	# move to install directory
	if [ ${MERLIN_INSTALL_DIR}/libmerlin_pg.a != ${MERLIN_BUILD_DIR}/lib/libmerlin_pg.a ] ; then \
		mv ${MERLIN_BUILD_DIR}/lib/libmerlin_pg.a ${MERLIN_INSTALL_DIR}/libmerlin_pg.a ; \
	fi


