/*

This source file exists in order to handle all merlin related mpi operations
and to keep as much as possible inside the merlin libary, instead of in individual
simulation files.

We work in SPMD land.

*/
#ifdef ENABLE_MPI

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
}


//Merlin worker node
//Merlin process support via choosing tags
void merlin_worker_node()
{
	bool do_work = true;
	MPI::Status merlin_mpi_status;
	MPI::Datatype MPI_Particle = merlin_mpi_create_particle();
	MPI::COMM_WORLD.Barrier();
	cout << "I am a node, waiting for data" << endl;
	while(do_work)
	{
		//Lets wait (blocking) and see what is going to be sent from the master node.
		cout << "Node waiting for probe" << endl;
		MPI::COMM_WORLD.Probe(0, MPI_ANY_TAG, merlin_mpi_status);
		
		tag = merlin_mpi_status.Get_tag();
		cout << "PROBE OK: " << tag << endl;

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
			cout << "NODE IN TRANSPORT" << endl;
			//First make a buffer to place the transport map inside
			//This is what will what will be applied.
			//RTMap *amap = new RTMap;
			
			//But also need a buffer to recv the one being sent, then convert.

			//Recv the map to apply
			//MPI::COMM_WORLD.Recv(transport_map, 1, MPI_Transport_Map, 0, 1, merlin_mpi_status);

			//Recv particles (unknown number)
			//We probe to see how many particles we will be getting - this blocks till a send occurs
			//only listens to sends from the master specifying transport particles.
			//MPI::COMM_WORLD.Probe(0, 1, merlin_mpi_status);

			//Find out the number of particles we are being given
			int particle_count = merlin_mpi_status.Get_count(MPI_Particle);
			cout << particle_count << endl;
			//We then make a suitable buffer for the particles
			PSvectorArray *particle_buffer_in = new PSvectorArray[particle_count];
			PSvectorArray *particle_buffer_out = new PSvectorArray[particle_count];
			cout << "good" << endl;
			//Recv + store particles
			//MPI::COMM_WORLD.Recv(buffer, int count, MPI_Particle, int source, int tag, MPI::Status);
			MPI::COMM_WORLD.Recv(&particle_buffer_in, particle_count, MPI_Particle, 0, 1, merlin_mpi_status);
			cout << "good" << endl;
			//We now have the map, and the particles required, so call the transport integrator function for each particle
			for (size_t n = 0; n<particle_count; n++)
			{
				//typedef std::vector<PSvector> PSvectorArray;
				//for_each(bunch.begin(),bunch.end(),ApplyMap(amap));
				cout << particle_buffer_in->size() << endl;
				//amap->Apply(particle_buffer_in);
				//cout << particle_buffer_in[n].v[0] << endl;
//cout << particle_buffer_in[n][0] << "\t" << particle_buffer_in[n][1] << "\t" << particle_buffer_in[n][2] << "\t" << particle_buffer_in[n][3] << "\t" << particle_buffer_in[n][4] << "\t" << particle_buffer_in[n][5] << endl;
			}
			//re-write array
			particle_buffer_out = particle_buffer_in;

			//Send back particles to the master node.
			MPI::COMM_WORLD.Send(&particle_buffer_out[0], particle_count, MPI_Particle, 0, 1);
			
			//Check all is ok?
			//end of this func -> (re-loop)
			//Clean up
			delete particle_buffer_in;
			delete particle_buffer_out;


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
/*
//Particle transport conversion
void merlin_mpi_get_particles(ParticleBunch& bunch)
{
	//PSvectorArray& ParticleBunch::GetParticles ()
	PSvectorArray tempbunch = bunch.GetParticles()
	
	//We can now access each particle with tempbunch[i]
	//This should be a PSvector, aka a particle.
	//individual coordinates should be tempbunch[i][0] for x
	//tempbunch[i][1] for x' etc
	//Since MPI_Particle is defined as 6 doubles, we can do the following:
	//MPI::COMM_WORLD.Send(tempbunch[i][0], nparticles, MPI_Particle, node, tag);
}
*/
MPI::Datatype merlin_mpi_create_particle()
{
	//The particle type MPI::Datatype
	
	MPI::Datatype MPI_Particle = MPI::DOUBLE.Create_contiguous(6);
	MPI_Particle.Commit();
	return MPI_Particle;
}

#endif
