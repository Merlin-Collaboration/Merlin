#ifndef _MERLIN_MPI_HPP_
#define _MERLIN_MPI_HPP_
#ifdef ENABLE_MPI
#include <mpi.h>

extern int size;
extern int rank;
extern int tag;
extern MPI::Datatype MPI_Particle;
extern MPI::Datatype MPI_Transport_Map;

void merlin_mpi_init(int argc, char* argv[]);

void merlin_worker_node();

//Create a particle type for network sends - 6 doubles
// x x' y y' ct dp
//See BeamModel/PSvector.h


#endif
#endif
