/*

This source file exists in order to handle all merlin related mpi operations
and to keep as much as possible inside the merlin libary, instead of in individual
simulation files.

We work in SPMD land.

*/

#include "merlin_mpi.hpp"
#include "BeamModel/PSvector.h"

using namespace std;
int size, rank, tag;

void merlin_mpi_init(int argc, char* argv[])
{
	//Start up MPI
	MPI::Init(argc, argv);

	//Find the total number of processes
	size = MPI::COMM_WORLD.Get_size();

	if (size <= 2)
	{
		cout << "At least 2 nodes are required currently." << endl;
		abort();
	}

	//find this processes rank
	rank = MPI::COMM_WORLD.Get_rank();
	
	//Also require some derived data types
	//The particle type	
	MPI::Datatype MPI_Particle = MPI::DOUBLE.Create_contiguous(6);
	MPI_Particle.Commit();

	//The transport map type	
	MPI::Datatype MPI_Transport_Map = MPI::DOUBLE.Create_contiguous(6);
	MPI_Transport_Map.Commit();
}


//Merlin worker node
//Merlin process support via choosing tags
void merlin_worker_node()
{
	bool do_work = true;
	MPI::Status merlin_mpi_status;
	cout << "I am a node, waiting for data" << endl;
	while(do_work)
	{
		//Lets wait (blocking) and see what is going to be sent from the master node.
		MPI::COMM_WORLD.Probe(0, MPI_ANY_TAG, merlin_mpi_status);
		tag = merlin_mpi_status.Get_tag();
		/*
		TAGS - these are entirely arbitrary choices for each merlin process.
		
		0	EXIT
		1	TRANSPORT
		2	COLLIMATION
		3	WAKEFIELDS
		
		*/

		if (tag == 0)		//EXIT
		{
			do_work = false;
			cout << "Node " << rank << " signaled to exit" << endl;
			break;
		}

		else if(tag == 1)	//TRANSPORT
		{
			///BeamDynamics/ParticleTracking/Integrators/TransportIntegrators.cpp

			//First make a buffer to place the transport map inside
			//This is what will what will be applied.
			RTMap *amap = new RTMap;
			
			//But also need a buffer to recv the one being sent, then convert.

			//Recv the map to apply
			MPI::COMM_WORLD.Recv(transport_map, 1, MPI_Transport_Map, 0, 1, merlin_mpi_status);

			//Recv particles (unknown number)
			//We probe to see how many particles we will be getting - this blocks till a send occurs
			//only listens to sends from the master specifying transport particles.
			MPI::COMM_WORLD.Probe(0, 1, merlin_mpi_status);

			//Find out the number of particles we are being given
			int particle_count = merlin_mpi_status.Get_count(MPI_Particle);

			//We then make a suitable buffer for the particles
			PSvector *particle_buffer_in = new PSvector[particle_count];
			PSvector *particle_buffer_out = new PSvector[particle_count];

			//Recv + store particles
			//MPI::COMM_WORLD.Recv(buffer, int count, MPI_Particle, int source, int tag, MPI::Status);
			MPI::COMM_WORLD.Recv(particle_buffer_in, particle_count, MPI_Particle, 0, 1, merlin_mpi_status);

			//We now have the map, and the particles required, so call the transport integrator function for each particle
			for (size_t n = 1; n<particle_buffer_in.size(); n++)
			{
				amap->Apply(particle_buffer_in);
			}
			//re-write array
			particle_buffer_out = particle_buffer_in;

			//Send back particles to the master node.
			MPI::COMM_WORLD.Send(particle_buffer_out, particle_count, MPI_Particle, 0, 1);
			
			//Check all is ok?
			//end of this func -> (re-loop)
			//Clean up
			delete particle_buffer_in;


			//Since collimation will be dropping particles, one can assume the order
			//in which particles are sent back to the master thread does not matter.
			//That is, for example the 10th particle in the initial array may not be
			//the 10th particle in the final array after tracking.
		}
		
		else if(tag == 2)	//COLLIMATION
		{
			//do_collimation();
		}
		
		else if(tag == 3)	//WAKEFIELDS
		{
			//do_wakefields();
		}
		else
		{
			//Should not get here, unknown tag/something went wrong
			//Some sort of error handling should exist, re-request?
			//For now, just break out
			cout << "Unknown MPI tag: " << tag << "\t Exiting node" << endl;
			do_work = false;
		}
	}

}

void merlin_mpi_finalize()
{
	//Always remember to clean up
	MPI::Finalize();
}
