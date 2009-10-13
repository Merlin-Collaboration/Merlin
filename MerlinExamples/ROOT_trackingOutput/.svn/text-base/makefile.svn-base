# Makefile for the new MERLIN Examples
# -------------------------------------
# 25.11.2008
# Dirk Kruecker, DESY
#
# The env variable ROOTSYS must be defined!
#
#

ifeq ($(strip $(ROOTSYS)),) 
$(error Please, define ROOTSYS!)
endif

CPP= g++

MERLIN_ROOT_DIR=${PWD}/../../Merlin
MERLIN_LIB_PATH=${PWD}/../../merlin-unix/lib
INCPATH= -I${MERLIN_ROOT_DIR} -I${ROOTSYS}/include
CCFLAGS2= -w -O3 ${INCPATH}
LDFLAGS2= `root-config --cflags --libs` -L${MERLIN_LIB_PATH}

MerlinROOT : MerlinROOT_main.o TeslaWakePotential.o TrackingOutputROOT.o 
	${CPP} ${LDFLAGS2} -o $@ $^ -lmerlin -lm

%.o : %.cpp
	${CPP} ${CCFLAGS2} -c  -o $@ $<
