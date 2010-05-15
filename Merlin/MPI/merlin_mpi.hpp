#ifndef _MERLIN_MPI_HPP_
#define _MERLIN_MPI_HPP_
#ifdef ENABLE_MPI
#include <mpi.h>

//int size,rank,tag;

//extern MPI::Datatype MPI_Transport_Map;

void merlin_mpi_init(int argc, char* argv[]);
MPI::Datatype merlin_mpi_create_particle();
void merlin_worker_node();
void merlin_mpi_finalize();

//Create a particle type for network sends - 6 doubles
// x x' y y' ct dp
//See BeamModel/PSvector.h

//Define the MPI_Particle type as a structure
//Should be easy to add extra fields to this, for macrocharge etc.
/*const int coords = 6;
int count = 1;
int blockcounts[1];
MPI::Aint offsets[1];
MPI::Datatype oldtypes[1];
offsets[0] = 0;
oldtypes[0] = MPI_DOUBLE;
blockcounts[0] = coords;

MPI::Datatype MPI_Particle = MPI::DOUBLE.Create_struct(count,blockcounts,offsets,oldtypes);
MPI_Particle.Commit();
*/
#endif
#endif
