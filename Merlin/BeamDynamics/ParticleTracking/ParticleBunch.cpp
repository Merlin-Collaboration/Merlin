/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:53 $
// $Revision: 1.7 $
// 
/////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <fstream>
#include <list>
#include <iterator>
// Transform3D
#include "EuclideanGeometry/Transform3D.h"
// PSvectorTransform3D
#include "BasicTransport/PSvectorTransform3D.h"
// ParticleBunch
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

using namespace std;

// needed by sort algorithm
inline bool operator<(const PSvector& p1, const PSvector& p2)
{
    return p1.ct()<p2.ct();
}

namespace {

using namespace ParticleTracking;

template<class V>
struct APSV1 {
    V& p;
    const double w;

    APSV1(V& p0, double n) : p(p0),w(1/n) {}
    void operator()(const PSvector& p1) {
        for(int i=0; i<6; i++) p[i]+=w*p1[i];
    }
};

struct APSV2 {
    const int u,v;
    Point2D X;
    const double w;

    APSV2(int i,int j,double n) : u(i),v(j),X(0,0),w(1/n){}
    void operator()(const PSvector& p) {
        X.x+=w*p[u];
        X.y+=w*p[v];
    }
};

struct PSVVAR1 {
    PSmoments& S;
    const double w;
    PSVVAR1(PSmoments& sig,double n):S(sig),w(1/n) {}
    void operator()(const PSvector& p) {
        for(int i=0; i<6; i++)
            for(int j=0; j<=i; j++)
                S(i,j)+=w*(p[i]-S[i])*(p[j]-S[j]);
    }
};

struct PSVVAR2 {
    PSmoments2D& S;
    const int u,v;
    const double w;

    PSVVAR2(int u1, int v1, double n, PSmoments2D& sig)
            :S(sig),u(u1),v(v1),w(1/n) {}

    void operator()(const PSvector& p) {
        double x=p[u]-S[0];
        double y=p[v]-S[1];
        S(0,0)+=w*x*x;
        S(0,1)+=w*x*y;
        S(1,1)+=w*y*y;
    }
};

template<class IT>
double Mean(IT F, IT L,int i)
{
    double s=0;
    double n=0;
    while(F!=L) {s+=(*F++)[i];n++;}
    return s/n;
}

template<class T>
inline void SortArray(std::vector<T>& array)
{
    sort(array.begin(),array.end());
}

template<class T>
inline void SortArray(std::list<T>& array)
{
    array.sort();
}

};

namespace ParticleTracking {

ParticleBunch::ParticleBunch (double P0, double Q, PSvectorArray& particles)
        : Bunch(P0,Q),qPerMP(Q/particles.size()),pArray(),init(false),coords(6)
{
    pArray.swap(particles);
}

ParticleBunch::ParticleBunch (double P0, double Q, std::istream& is)
        : Bunch(P0,Q),init(false),coords(6)
{
    PSvector p;
    while(is>>p)
        push_back(p);

    qPerMP = Q/size();
}

ParticleBunch::ParticleBunch (double P0, double Qm)
        : Bunch(P0,Qm),qPerMP(Qm),init(false),coords(6)
{}

double ParticleBunch::GetTotalCharge () const
{
    return qPerMP*size();
}

PSmoments& ParticleBunch::GetMoments (PSmoments& sigma) const
{
    sigma.zero();
    for_each(begin(),end(),APSV1<PSmoments>(sigma,size()));
    for_each(begin(),end(),PSVVAR1(sigma,size()));
    return sigma;
}

PSmoments2D& ParticleBunch::GetProjectedMoments (PScoord u, PScoord v, PSmoments2D& sigma) const
{
    Point2D X = for_each(begin(),end(),APSV2(u,v,size())).X;
    sigma[0]=X.x;
    sigma[1]=X.y;
    for_each(begin(),end(),PSVVAR2(u,v,size(),sigma));
    return sigma;
}

PSvector& ParticleBunch::GetCentroid (PSvector& p) const
{
    p.zero();
    for_each(begin(),end(),APSV1<PSvector>(p,size()));
    return p;
}

std::pair<double,double> ParticleBunch::GetMoments(PScoord i) const
{
    double u = Mean(begin(),end(),i);
    double v=0;
    for(const_iterator p = begin(); p!=end(); p++) {
        double x = (*p)[i]-u;
        v+=x*x;
    }

    return make_pair(u,sqrt(v/size()));
}

Point2D ParticleBunch::GetProjectedCentroid (PScoord u, PScoord v) const
{
    return for_each(begin(),end(),APSV2(u,v,size())).X;
}

double ParticleBunch::AdjustRefMomentumToMean ()
{
    return AdjustRefMomentum( Mean(begin(),end(),ps_DP) );
}

double ParticleBunch::AdjustRefMomentum (double dpp)
{
    double onePlusDpp = 1+dpp;

    for(iterator p=begin(); p!=end(); p++)
        p->dp() = (p->dp()-dpp)/onePlusDpp;

    double P0 = onePlusDpp*GetReferenceMomentum();
    SetReferenceMomentum(P0);
    return P0;
}

double ParticleBunch::AdjustRefTimeToMean ()
{
    double meanct = Mean(begin(),end(),ps_CT);
    for(iterator p=begin(); p!=end(); p++)
        (*p).ct()-=meanct;

    double CT = GetReferenceTime()-meanct;
    SetReferenceTime(CT);
    return CT;
}

Histogram& ParticleBunch::ProjectDistribution (PScoord axis, Histogram& hist) const
{
    // TODO:
    return hist;
}


bool ParticleBunch::ApplyTransformation (const Transform3D& t)
{
    if(!t.isIdentity())
        PSvectorTransform3D(t).Apply(pArray);
    return true;
}

void ParticleBunch::SortByCT ()
{
    //	pArray.sort();
    SortArray(pArray);
}

void ParticleBunch::Output (std::ostream& os) const
{
    //	std::copy(begin(),end(),ostream_iterator<PSvector>(os));
//    int oldp=os.precision(10);
    int oldp=os.precision(16);
    ios_base::fmtflags oflg = os.setf(ios::scientific,ios::floatfield);
    for(PSvectorArray::const_iterator p = begin(); p!=end(); p++) {
//        os<<std::setw(24)<<GetReferenceTime();
//        os<<std::setw(24)<<GetReferenceMomentum();
        os<<std::setw(35)<<GetReferenceTime();
        os<<std::setw(35)<<GetReferenceMomentum();
        for(size_t k=0; k<6; k++)
            os<<std::setw(35)<<(*p)[k];//20
        os<<endl;
    }
    os.precision(oldp);
    os.flags(oflg);
}

void ParticleBunch::Input (double Q, std::istream& is)
{
    double reftime, refmom;
    PSvector p;
    while(is>>reftime>>refmom>>p){
        push_back(p);

 	   SetReferenceTime(reftime);
    	SetReferenceMomentum(refmom);
}
    qPerMP = Q/size();
}



void ParticleBunch::SetCentroid (const Particle& x0)
{
    PSvector x;
    GetCentroid(x);
    x-=x0;
    for(PSvectorArray::iterator p = begin(); p!=end(); p++)
        *p-=x;
}

}; // end namespace ParticleTracking

//MPI code
#ifdef ENABLE_MPI

void ParticleBunch::MPI_Initialize()
{
	//Check of the MPI runtime has started
	if (!MPI::Is_initialized())
	{
		//If not, start it.
		MPI::Init();
	}

	/*
	We want MPI to throw C++ exceptions that we can then deal with. For example:
	
		try
		{
			MPI::something();
		}
		catch(MPI::Exception fail)
		{
			cout << fail.Get_error_string() << endl;
			cout << fail.Get_error_code() << endl;
			MPI::COMM_WORLD.Abort();
		}

	*/
	MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
	
	//Total number of processors in the cluster
	MPI_size = MPI::COMM_WORLD.Get_size();

	//Make sure we have 1 master, 2 nodes until the particle distribution is "fixed"
	if (MPI_size <= 2)
	{
		cout << "At least 2 nodes are currently required." << endl;
		MPI::COMM_WORLD.Abort(1);
	}

	//find this processes rank
	MPI_rank = MPI::COMM_WORLD.Get_rank();

	//Create the particle type
	Create_MPI_particle();
}

void ParticleBunch::Create_MPI_particle()
{
	//The particle type MPI::Datatype
	int count = 1;
	int blockcounts[1];
	MPI::Aint offsets[1];
	MPI::Datatype oldtypes[1];
	offsets[0] = 0;
	oldtypes[0] = MPI_DOUBLE;
	blockcounts[0] = coords;

	//We now create the new datatype
	try
	{
		MPI_Particle = MPI::DOUBLE.Create_struct(count,blockcounts,offsets,oldtypes);
	}
	catch(MPI::Exception fail)
	{
		cout << "MPI Particle type creation fail in Create_MPI_particle() - Create_struct()" << endl;
		cout << "Error :" << fail.Get_error_string() << "\tNumber: " << fail.Get_error_code() << endl;
		MPI::COMM_WORLD.Abort(1);
	}

	//And commit it.
	try
	{
		MPI_Particle.Commit();
	}
	catch(MPI::Exception fail)
	{
		cout << "MPI Particle type creation fail in Create_MPI_particle() - Commit()" << endl;
		cout << "Error :" << fail.Get_error_string() << "\tNumber: " << fail.Get_error_code() << endl;
		MPI::COMM_WORLD.Abort(1);
	}

}

void ParticleBunch::MPI_Finalize()
{
	//Clean up at exit.
	//Free the created Particle type
	MPI_Particle.Free();
	//And finalize the MPI process
	MPI::Finalize();
}

void ParticleBunch::master_recv_particles_from_nodes()
{
	PSvector Particle_buffer;
	for(int n = 1; n<MPI_size; n++)
	{
               	//We probe the incomming send to see how many particles we are recieving.
                MPI::COMM_WORLD.Probe(n, MPI_ANY_TAG, MPI_status);

                //We find out from the probe the number of particles we have been sent
                int recv_count = MPI_status.Get_count(MPI_Particle);

                //We make a suitable buffer
		try
		{
			particle_recv_buffer = new double[recv_count*coords];
		}
		catch(std::bad_alloc)
		{
			cout << "bad_alloc in master_recv_particles_from_nodes() on node " << MPI_rank << " - recv count: " << recv_count << "\t coords: " << coords << endl;
			MPI::COMM_WORLD.Abort(1);
		}

       	        //We now do the Recv for real
                MPI::COMM_WORLD.Recv(particle_recv_buffer, recv_count, MPI_Particle, n, 1, MPI_status);

       	        //Put the recv buffer into the particle array
               	for (int i = 0; i< recv_count; i++)
                {
       	                Particle_buffer[0] = particle_recv_buffer[(i*coords)+0];
                        Particle_buffer[1] = particle_recv_buffer[(i*coords)+1];
                        Particle_buffer[2] = particle_recv_buffer[(i*coords)+2];
                        Particle_buffer[3] = particle_recv_buffer[(i*coords)+3];
                        Particle_buffer[4] = particle_recv_buffer[(i*coords)+4];
                        Particle_buffer[5] = particle_recv_buffer[(i*coords)+5];

                        //Push back each particle onto the current master node bunch.
                        push_back(Particle_buffer);
       	        }

	delete [] particle_recv_buffer;
	particle_recv_buffer = NULL;
	}//Merge bunches now finished.
}


void ParticleBunch::node_recv_particles_from_master()
{
	//We probe the incomming send to see how many particles we are recieving.
	MPI::COMM_WORLD.Probe(0, MPI_ANY_TAG, MPI_status);
	
	//We find out from the probe the number of particles we have been sent
	int recv_count = MPI_status.Get_count(MPI_Particle);

	//We make a suitable buffer
	try
	{
		particle_recv_buffer = new double[recv_count*coords];
	}
	catch(std::bad_alloc)
	{
		cout << "bad_alloc in node_recv_particles_from_master() on node " << MPI_rank << " - recv count: " << recv_count << "\t coords: " << coords << endl;
		MPI::COMM_WORLD.Abort(1);
	}

	//We now do the Recv for real, 
	MPI::COMM_WORLD.Recv(particle_recv_buffer, recv_count, MPI_Particle, 0, 1, MPI_status);

	//We now have our particle data stored in an array
	//Must now convert them into particles, then a bunch.
	PSvector Particle_buffer;

	//Put the recv buffer into the particle array
	for (int i = 0; i< recv_count; i++)
	{
		Particle_buffer[0] = particle_recv_buffer[(i*coords)+0];
		Particle_buffer[1] = particle_recv_buffer[(i*coords)+1];
		Particle_buffer[2] = particle_recv_buffer[(i*coords)+2];
		Particle_buffer[3] = particle_recv_buffer[(i*coords)+3];
		Particle_buffer[4] = particle_recv_buffer[(i*coords)+4];
		Particle_buffer[5] = particle_recv_buffer[(i*coords)+5];

		//Push back each particle onto the cleared currentBunch
		push_back(Particle_buffer);
	}

	delete [] particle_recv_buffer;
	particle_recv_buffer = NULL;
}

void ParticleBunch::node_send_particles_to_master()
{
	//Here we are a "worker" node, and will send particles to the master for the collective calculation
	//send bunches
	int particle_count = size();
	
	//We make an array to put the particle bunch data into.
	//It needs to be continuous in memory to work nicely with MPI
	//This is crude, but works.
	try
	{
		particle_send_buffer = new double[particle_count*coords];
	}
	catch(std::bad_alloc)
	{
		cout << "bad_alloc in node_send_particles_to_master() on node " << MPI_rank << " - particle count: " << particle_count << "\t coords: " << coords << endl;
		MPI::COMM_WORLD.Abort(1);
	}

	for (int n=0; n<particle_count; n++)
	{
		particle_send_buffer[(n*coords)+0] = GetParticles()[n].x();
		particle_send_buffer[(n*coords)+1] = GetParticles()[n].xp();
		particle_send_buffer[(n*coords)+2] = GetParticles()[n].y();
		particle_send_buffer[(n*coords)+3] = GetParticles()[n].yp();
		particle_send_buffer[(n*coords)+4] = GetParticles()[n].ct();
		particle_send_buffer[(n*coords)+5] = GetParticles()[n].dp();
	}
		
	//Send everything to the master node
	MPI::COMM_WORLD.Send(&particle_send_buffer[0], particle_count, MPI_Particle, 0, 1);

	//Whilst the master node is busy calculating the Wakefield, we may as well use this time to clear the currentBunch.
	//This will be refilled with new particles.
	clear();
	
	delete [] particle_send_buffer;
	particle_send_buffer = NULL;
}

void ParticleBunch::master_send_particles_to_nodes()
{
	//Now save the bunch size
	int particle_count = size();

	//The bunch must now be convered into a format suitable for sending again.
	//The leftover extra particles can be added to the "local" machine bunch.
	int particles_per_node = particle_count/MPI_size;
	int remaining_particles = particle_count % MPI_size;

	try
	{
		particle_send_buffer = new double[particle_count*coords];
	}
	catch(std::bad_alloc)
	{
		cout << "bad_alloc in master_send_particles_to_nodes() on node " << MPI_rank << "  - particle count: " << particle_count << "\t coords: " << coords << endl;
		MPI::COMM_WORLD.Abort(1);
	}

	for (int n=0; n<particle_count; n++)
	{
		particle_send_buffer[(n*coords)+0] = GetParticles()[n].x();
		particle_send_buffer[(n*coords)+1] = GetParticles()[n].xp();
		particle_send_buffer[(n*coords)+2] = GetParticles()[n].y();
		particle_send_buffer[(n*coords)+3] = GetParticles()[n].yp();
		particle_send_buffer[(n*coords)+4] = GetParticles()[n].ct();
		particle_send_buffer[(n*coords)+5] = GetParticles()[n].dp();
	}
	//We now have the full bunch in a particle buffer.
	//This must be sliced, and particles sent to each node

	//Resend
	// 1 -> total nodes
	for(int n = 0; n<(MPI_size-1); n++)
	{
		MPI::COMM_WORLD.Send(&particle_send_buffer[(particles_per_node*n*coords)], particles_per_node, MPI_Particle, (n+1), 1);
	}

	//Finally the local bunch must be reconstructed from the remaining particles.
	//First it must be cleared
	clear();

	PSvector Particle_buffer;
	//And then refilled.
	for (int i = (particles_per_node*(MPI_size-1)); i<particle_count; i++)
	{
		Particle_buffer[0] = particle_send_buffer[(i*coords)+0];
		Particle_buffer[1] = particle_send_buffer[(i*coords)+1];
		Particle_buffer[2] = particle_send_buffer[(i*coords)+2];
		Particle_buffer[3] = particle_send_buffer[(i*coords)+3];
		Particle_buffer[4] = particle_send_buffer[(i*coords)+4];
		Particle_buffer[5] = particle_send_buffer[(i*coords)+5];

		//Push back each particle onto the cleared currentBunch
		push_back(Particle_buffer);
	}
	delete [] particle_send_buffer;
	particle_send_buffer = NULL;

}

//Gather particle function: All particles on the nodes are moved to the bunch on the master node.
void ParticleBunch::gather()
{

if (init == false)
{
	MPI_Initialize();
	init = true;
}
	MPI::COMM_WORLD.Barrier();

	if (MPI_rank == 0)
	{
		master_recv_particles_from_nodes();
	}

	else
	{
		node_send_particles_to_master();
	}
}

//Particle distribution function: Here all particles on the master are distributed between the nodes/
void ParticleBunch::distribute()
{

if (init == false)
{
	MPI_Initialize();
	init = true;
}
	MPI::COMM_WORLD.Barrier();
	if (MPI_rank == 0)
	{
		master_send_particles_to_nodes();
	}

	else
	{
		node_recv_particles_from_master();
	}
}

//Destructor: Only enabled with MPI - Will clean up and finalize.
ParticleBunch::~ParticleBunch ()
{
	MPI_Finalize();
}
#endif	//end MPI code
