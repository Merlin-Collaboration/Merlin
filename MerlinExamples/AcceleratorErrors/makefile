# Makefile for the new MERLIN Examples
# -------------------------------------
# 25.11.2008
# Dirk Kruecker, DESY
#

CPP= g++
MERLIN_ROOT_DIR=${PWD}/../../Merlin
MERLIN_LIB_PATH=${PWD}/../../merlin-unix/lib
INCPATH= -I${MERLIN_ROOT_DIR}
CCFLAGS= -w -O3 ${INCPATH}
LDFLAGS= -L${MERLIN_LIB_PATH} -static 

AccErr : main.o
	${CPP} ${LDFLAGS} -o $@ $^ -lmerlin -lm

%.o : %.cpp
	${CPP} ${CCFLAGS} -c  -o $@ $<
