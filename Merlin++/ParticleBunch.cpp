/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <algorithm>
#include <fstream>
#include <list>
#include <iterator>
#include "Transform3D.h"
#include "PSvectorTransform3D.h"
#include "ParticleBunch.h"
#include "NormalTransform.h"
#include "MatrixMaps.h"
#include "ParticleDistributionGenerator.h"
#include "BeamData.h"
#include "BunchFilter.h"

#ifdef MERLIN_PROFILE
#include "MerlinProfile.h"
#endif

using namespace std;

// needed by sort algorithm
inline bool operator<(const PSvector& p1, const PSvector& p2)
{
	return p1.ct() < p2.ct();
}

namespace
{

using namespace ParticleTracking;

template<class V>
struct APSV1
{
	V& p;
	double w;

	APSV1(V& p0, double n) :
		p(p0), w(0)
	{
	}
	void operator()(const PSvector& p1)
	{
		w++;
		for(int i = 0; i < 6; i++)
		{
			p[i] += (p1[i] - p[i]) / w;
		}
	}

};

struct APSV2
{
	const int u, v;
	Point2D X;
	const double w;

	APSV2(int i, int j, double n) :
		u(i), v(j), X(0, 0), w(1 / n)
	{
	}
	void operator()(const PSvector& p)
	{
		X.x += w * p[u];
		X.y += w * p[v];
	}

};

struct PSVVAR1
{
	PSmoments& S;
	double w;
	PSVVAR1(PSmoments& sig, double n) :
		S(sig), w(0)
	{
	}
	void operator()(const PSvector& p)
	{
		w++;
		for(int i = 0; i < 6; i++)
			for(int j = 0; j <= i; j++)
			{
				S(i, j) += ((p[i] - S[i]) * (p[j] - S[j]) - S(i, j)) / w;
			}
	}

};

struct PSVVAR2
{
	PSmoments2D& S;
	const int u, v;
	const double w;

	PSVVAR2(int u1, int v1, double n, PSmoments2D& sig) :
		S(sig), u(u1), v(v1), w(1 / n)
	{
	}

	void operator()(const PSvector& p)
	{
		double x = p[u] - S[0];
		double y = p[v] - S[1];
		S(0, 0) += w * x * x;
		S(0, 1) += w * x * y;
		S(1, 1) += w * y * y;
	}

};

template<class IT>
double Mean(IT F, IT L, int i)
{
	double s = 0;
	double n = 0;
	while(F != L)
	{
		s += ((*F++)[i] - s) / (++n);
	}
	return s;
}

template<class T>
inline void SortArray(std::vector<T>& array)
{
	sort(array.begin(), array.end());
}

template<class T>
inline void SortArray(std::list<T>& array)
{
	array.sort();
}

} //end namespace

namespace ParticleTracking
{

ParticleBunch::ParticleBunch(double P0, double Q, PSvectorArray& particles) :
	Bunch(P0, Q), init(false), coords((int) sizeof(PSvector) / sizeof(double)), ScatteringPhysicsModel(0), qPerMP(Q
		/ particles.size()), pArray()
{
	pArray.swap(particles);
}

ParticleBunch::ParticleBunch(double P0, double Q, std::istream& is) :
	Bunch(P0, Q), init(false), coords((int) sizeof(PSvector) / sizeof(double)), ScatteringPhysicsModel(0)
{
	PSvector p;
	while(is >> p)
	{
		push_back(p);
	}

	qPerMP = Q / size();
}

ParticleBunch::ParticleBunch(double P0, double Qm) :
	Bunch(P0, Qm), init(false), coords((int) sizeof(PSvector) / sizeof(double)), ScatteringPhysicsModel(0), qPerMP(Qm)
{
}

ParticleBunch::ParticleBunch(size_t np, const ParticleDistributionGenerator& generator, const BeamData& beam,
	ParticleBunchFilter* filter) :
	ParticleBunch(beam.p0, beam.charge)
{
	RMtrx M;
	M.R = NormalTransform(beam);

	// The first particle is *always* the centroid particle
	PSvector p;
	p.x() = beam.x0;
	p.xp() = beam.xp0;
	p.y() = beam.y0;
	p.yp() = beam.yp0;
	p.dp() = 0;
	p.ct() = beam.ct0;
	p.type() = -1.0;
	p.location() = -1.0;
	p.id() = 0;
	p.sd() = 0.0;
	pArray.push_back(p);

	size_t i = 1;
	size_t filtered = 0;
	while(i < np)
	{
		p = generator.GenerateFromDistribution();

		// apply emittance
		p.x() *= sqrt(beam.emit_x);
		p.xp() *= sqrt(beam.emit_x);
		p.y() *= sqrt(beam.emit_y);
		p.yp() *= sqrt(beam.emit_y);
		p.dp() *= sqrt(beam.sig_dp);
		p.ct() *= sqrt(beam.sig_z);

		// Apply Courant-Snyder
		M.Apply(p);

		p += pArray.front(); // add centroid

		p.type() = -1.0;
		p.location() = -1.0;
		p.id() = i;
		p.sd() = 0.0;

		if(filter == nullptr || filter->Apply(p))
		{
			pArray.push_back(p);
			i++;
		}
		else
		{
			filtered += 1;
		}

		if(filtered == 10000 && i == 1)
		{
			cout << "WARNING: In ParticleBunch::ParticleBunch() ParticleBunchFilter has skipped over " << filtered
				 << " particles without allowing any." << endl;
		}
	}
	qPerMP = beam.charge / size();
}

double ParticleBunch::GetTotalCharge() const
{
	return qPerMP * size();
}

PSmoments& ParticleBunch::GetMoments(PSmoments& sigma) const
{
	sigma.zero();
	for_each(begin(), end(), APSV1<PSmoments>(sigma, size()));
	for_each(begin(), end(), PSVVAR1(sigma, size()));
	return sigma;
}

PSmoments2D& ParticleBunch::GetProjectedMoments(PScoord u, PScoord v, PSmoments2D& sigma) const
{
	Point2D X = for_each(begin(), end(), APSV2(u, v, size())).X;
	sigma[0] = X.x;
	sigma[1] = X.y;
	for_each(begin(), end(), PSVVAR2(u, v, size(), sigma));
	return sigma;
}

PSvector& ParticleBunch::GetCentroid(PSvector& p) const
{
	p.zero();
	for_each(begin(), end(), APSV1<PSvector>(p, size()));
	return p;
}

std::pair<double, double> ParticleBunch::GetMoments(PScoord i) const
{
	double u = Mean(begin(), end(), i);
	double v = 0;
	double w = 0;
	for(const_iterator p = begin(); p != end(); p++)
	{
		double x = (*p)[i] - u;
		v += (x * x - v) / (++w);
	}

	return make_pair(u, sqrt(v));
}

Point2D ParticleBunch::GetProjectedCentroid(PScoord u, PScoord v) const
{
	return for_each(begin(), end(), APSV2(u, v, size())).X;
}

double ParticleBunch::AdjustRefMomentumToMean()
{
	return AdjustRefMomentum(Mean(begin(), end(), ps_DP));
}

double ParticleBunch::AdjustRefMomentum(double dpp)
{
	double onePlusDpp = 1 + dpp;

	for(iterator p = begin(); p != end(); p++)
	{
		p->dp() = (p->dp() - dpp) / onePlusDpp;
	}

	double P0 = onePlusDpp * GetReferenceMomentum();
	SetReferenceMomentum(P0);
	return P0;
}

double ParticleBunch::AdjustRefTimeToMean()
{
	double meanct = Mean(begin(), end(), ps_CT);
	for(iterator p = begin(); p != end(); p++)
	{
		(*p).ct() -= meanct;
	}

	double CT = GetReferenceTime() - meanct;
	SetReferenceTime(CT);
	return CT;
}

Histogram& ParticleBunch::ProjectDistribution(PScoord axis, Histogram& hist) const
{
	// TODO:
	return hist;
}

bool ParticleBunch::ApplyTransformation(const Transform3D& t)
{
	if(!t.isIdentity())
	{
		PSvectorTransform3D(t).Apply(pArray);
	}
	return true;
}

void ParticleBunch::SortByCT()
{
	SortArray(pArray);
}

void ParticleBunch::Output(std::ostream& os) const
{
	Output(os, true);
}

void ParticleBunch::Output(std::ostream& os, bool show_header) const
{
	int oldp = os.precision(16);
	ios_base::fmtflags oflg = os.setf(ios::scientific, ios::floatfield);
	if(show_header)
	{
		os << "#T P0 X XP Y YP CT DP" << std::endl;
	}
	for(PSvectorArray::const_iterator p = begin(); p != end(); p++)
	{
		os << std::setw(35) << GetReferenceTime();
		os << std::setw(35) << GetReferenceMomentum();
		for(size_t k = 0; k < 6; k++)
		{
			os << std::setw(35) << (*p)[k];
		}
		os << endl;
	}
	os.precision(oldp);
	os.flags(oflg);
}

void ParticleBunch::OutputIndexParticle(std::ostream& os, int index) const
{
	int oldp = os.precision(16);
	ios_base::fmtflags oflg = os.setf(ios::scientific, ios::floatfield);

	PSvectorArray::const_iterator p = begin() + index;

	for(size_t k = 0; k < 6; k++)
	{
		os << std::setw(35) << (*p)[k];
	}
	os << endl;
	os.precision(oldp);
	os.flags(oflg);
}

void ParticleBunch::Input(double Q, std::istream& is)
{
	double reftime, refmom;
	PSvector p;
	string line;
	while(getline(is, line))
	{
		if(line == "" || line[0] == '#')
		{
			continue;
		}
		istringstream liness(line);
		liness >> reftime >> refmom >> p;
		push_back(p);

		SetReferenceTime(reftime);
		SetReferenceMomentum(refmom);
	}
	qPerMP = Q / size();
}

void ParticleBunch::SetCentroid()
{
	PSvector x;
	GetCentroid(x);
	x -= FirstParticle();
	for(PSvectorArray::iterator p = begin() + 1; p != end(); p++)
	{
		p->x() -= x.x();
		p->xp() -= x.xp();
		p->y() -= x.y();
		p->yp() -= x.yp();
		p->dp() -= x.dp();
		p->ct() -= x.ct();
	}
}

void ParticleBunch::SetCentroid(const Particle& x0)
{
	FirstParticle() = x0;
	SetCentroid();
}

bool ParticleBunch::IsStable() const
{
	return true;
}

double ParticleBunch::GetParticleMass() const
{
	return 0;
}

double ParticleBunch::GetParticleMassMeV() const
{
	return 0;
}

double ParticleBunch::GetParticleLifetime() const
{
	return 0;
}

//MPI code
#ifdef ENABLE_MPI

void ParticleBunch::MPI_Initialize()
{
	//Check of the MPI runtime has started
	if(!MPI::Is_initialized())
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
	/*
	    //Make sure we have 1 master, 2 nodes until the particle distribution is "fixed"
	    if (MPI_size <= 2)
	    {
	        cout << "At least 2 nodes are currently required." << endl;
	        MPI::COMM_WORLD.Abort(1);
	    }
	 */
	//find this processes rank
	MPI_rank = MPI::COMM_WORLD.Get_rank();

	//Find the processor/hostname
//	MPI::Get_processor_name(proc_name,char_array_len);

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
		MPI_Particle = MPI::DOUBLE.Create_struct(count, blockcounts, offsets, oldtypes);
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
	for(int n = 1; n < MPI_size; n++)
	{
		//We probe the incomming send to see how many particles we are recieving.
		MPI::COMM_WORLD.Probe(n, MPI_ANY_TAG, MPI_status);

		//We find out from the probe the number of particles we have been sent
		int recv_count = MPI_status.Get_count(MPI_Particle);

		//We make a suitable buffer
		try
		{
			particle_recv_buffer = new double[recv_count * coords];
		}
		catch(std::bad_alloc)
		{
			cout << "bad_alloc in master_recv_particles_from_nodes() on node " << MPI_rank << " - recv count: "
				 << recv_count << "\t coords: " << coords << endl;
			MPI::COMM_WORLD.Abort(1);
		}

		//We now do the Recv for real
		MPI::COMM_WORLD.Recv(particle_recv_buffer, recv_count, MPI_Particle, n, 1, MPI_status);

		//Put the recv buffer into the particle array
		for(int i = 0; i < recv_count; i++)
		{
			for(int j = 0; j < coords; j++)
			{
				Particle_buffer[j] = particle_recv_buffer[(i * coords) + j];
			}
			//Push back each particle onto the current master node bunch.
			push_back(Particle_buffer);
		}

		delete[] particle_recv_buffer;
		particle_recv_buffer = NULL;
	} //Merge bunches now finished.
}

void ParticleBunch::node_recv_particles_from_master()
{
	//We clear the currentBunch.
	//This will be refilled with new particles.
	clear();

	//We probe the incomming send to see how many particles we are recieving.
	MPI::COMM_WORLD.Probe(0, MPI_ANY_TAG, MPI_status);

	//We find out from the probe the number of particles we have been sent
	int recv_count = MPI_status.Get_count(MPI_Particle);

	//We make a suitable buffer
	try
	{
		particle_recv_buffer = new double[recv_count * coords];
	}
	catch(std::bad_alloc)
	{
		cout << "bad_alloc in node_recv_particles_from_master() on node " << MPI_rank << " - recv count: "
			 << recv_count << "\t coords: " << coords << endl;
		MPI::COMM_WORLD.Abort(1);
	}

	//We now do the Recv for real,
	MPI::COMM_WORLD.Recv(particle_recv_buffer, recv_count, MPI_Particle, 0, 1, MPI_status);

	//We now have our particle data stored in an array
	//Must now convert them into particles, then a bunch.
	PSvector Particle_buffer;

	//Put the recv buffer into the particle array
	for(int i = 0; i < recv_count; i++)
	{
		for(int j = 0; j < coords; j++)
		{
			Particle_buffer[j] = particle_recv_buffer[(i * coords) + j];
		}

		//Push back each particle onto the cleared currentBunch
		push_back(Particle_buffer);
	}

	delete[] particle_recv_buffer;
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
		particle_send_buffer = new double[particle_count * coords];
	}
	catch(std::bad_alloc)
	{
		cout << "bad_alloc in node_send_particles_to_master() on node " << MPI_rank << " - particle count: "
			 << particle_count << "\t coords: " << coords << endl;
		MPI::COMM_WORLD.Abort(1);
	}

	//Alternate method - will only work with a vector
	//memcpy(particle_send_buffer,&input[0],sizeof(double)*coords*size());

	int n = 0;
	for(PSvectorArray::iterator ip = begin(); ip != end(); ip++)
	{
		//for (int j=0; j<coords; j++)
		for(int j = 0; j < 6; j++)
		{
			particle_send_buffer[(n * coords) + 0] = ip->x();
			particle_send_buffer[(n * coords) + 1] = ip->xp();
			particle_send_buffer[(n * coords) + 2] = ip->y();
			particle_send_buffer[(n * coords) + 3] = ip->yp();
			particle_send_buffer[(n * coords) + 4] = ip->ct();
			particle_send_buffer[(n * coords) + 5] = ip->dp();
		}
		n++;
	}

	//Send everything to the master node
	MPI::COMM_WORLD.Send(&particle_send_buffer[0], particle_count, MPI_Particle, 0, 1);

	delete[] particle_send_buffer;
	particle_send_buffer = NULL;
}

void ParticleBunch::master_send_particles_to_nodes()
{
	//Now save the bunch size
	int particle_count = size();

	//The bunch must now be converted into a format suitable for sending again.
	//The leftover extra particles can be added to the "local" machine bunch.
	int particles_per_node = particle_count / MPI_size;
	int remaining_particles = particle_count % MPI_size;

	try
	{
		particle_send_buffer = new double[particle_count * coords];
	}
	catch(std::bad_alloc)
	{
		cout << "bad_alloc in master_send_particles_to_nodes() on node " << MPI_rank << "  - particle count: "
			 << particle_count << "\t coords: " << coords << endl;
		MPI::COMM_WORLD.Abort(1);
	}

	//Test new code
	int n = 0;
	for(PSvectorArray::iterator ip = begin(); ip != end(); ip++)
	{
		//for (int j=0; j<coords; j++)
		for(int j = 0; j < 6; j++)
		{
			particle_send_buffer[(n * coords) + 0] = ip->x();
			particle_send_buffer[(n * coords) + 1] = ip->xp();
			particle_send_buffer[(n * coords) + 2] = ip->y();
			particle_send_buffer[(n * coords) + 3] = ip->yp();
			particle_send_buffer[(n * coords) + 4] = ip->ct();
			particle_send_buffer[(n * coords) + 5] = ip->dp();
		}
		n++;
	}

	//We now have the full bunch in a particle buffer.
	//This must be sliced, and particles sent to each node

	//Resend
	// 1 -> total nodes
	for(int n = 0; n < (MPI_size - 1); n++)
	{
		MPI::COMM_WORLD.Send(&particle_send_buffer[(particles_per_node * n * coords)], particles_per_node, MPI_Particle,
			(n + 1), 1);
	}

	//Finally the local bunch must be reconstructed from the remaining particles.
	//First it must be cleared
	clear();

	PSvector Particle_buffer;
	//And then refilled.
	for(int i = (particles_per_node * (MPI_size - 1)); i < particle_count; i++)
	{
		for(int j = 0; j < coords; j++)
		{
			Particle_buffer[j] = particle_send_buffer[(i * coords) + j];
		}

		//Push back each particle onto the cleared currentBunch
		push_back(Particle_buffer);
	}
	delete[] particle_send_buffer;
	particle_send_buffer = NULL;

}

//Gather particle function: All particles on the nodes are moved to the bunch on the master node.
void ParticleBunch::gather()
{
//Really must switch this to MPI::COMM_WORLD.Gatherv()
//http://www.open-mpi.org/doc/v1.3/man3/MPI_Gatherv.3.php
//Will involve a call to size() and sending this to the master first from each node
//MPI::COMM_WORLD.Gather() on 1 int each.

#ifdef MERLIN_PROFILE
	MerlinProfile::AddProcess("GATHER");
	MerlinProfile::StartProcessTimer("GATHER");
#endif

	if(init == false)
	{
		MPI_Initialize();
		init = true;
	}
	MPI::COMM_WORLD.Barrier();
	if(MPI_rank == 0)
	{
		master_recv_particles_from_nodes();
	}

	else
	{
		node_send_particles_to_master();
	}

#ifdef MERLIN_PROFILE
	MerlinProfile::EndProcessTimer("GATHER");
#endif
}

//Particle distribution function: Here all particles on the master are distributed between the nodes/
void ParticleBunch::distribute()
{
#ifdef MERLIN_PROFILE
	MerlinProfile::AddProcess("SCATTER");
	MerlinProfile::StartProcessTimer("SCATTER");
#endif

//Should use MPI::COMM_WORLD.Scatterv() as per the gather
	Check_MPI_init();
	MPI::COMM_WORLD.Barrier();

	//Set up
	if(MPI_rank == 0)
	{
		master_send_particles_to_nodes();
	}
	else
	{
		node_recv_particles_from_master();
	}

#ifdef MERLIN_PROFILE
	MerlinProfile::EndProcessTimer("SCATTER");
#endif

}

void ParticleBunch::SendReferenceMomentum()
{
	Check_MPI_init();

	MPI::COMM_WORLD.Barrier();
	if(MPI_rank == 0)
	{
		//Get reference momentum
		AdjustRefMomentumToMean();
		double refP = GetReferenceMomentum();

		//Send to nodes
		for(int n = 1; n < MPI_size; n++)
		{
			MPI::COMM_WORLD.Send(&refP, 1, MPI::DOUBLE, n, 1);
		}
	}
	else
	{
		//Make a buffer
		double refP;

		//Recv new momentum
		MPI::COMM_WORLD.Recv(&refP, 1, MPI::DOUBLE, 0, 1, MPI_status);

		//Set new ref momentum
		SetReferenceMomentum(refP);
	}

}

//Destructor: Only enabled with MPI - Will clean up and finalize.
ParticleBunch::~ParticleBunch()
{
	//MPI_Finalize();
}

void ParticleBunch::Check_MPI_init()
{
	if(init == false)
	{
		MPI_Initialize();
		init = true;
	}
}
#endif  //end MPI code

} // end namespace ParticleTracking
