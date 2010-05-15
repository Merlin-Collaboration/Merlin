// Modified by D.Kruecker 18.2.2008
// to be used as base class for other wakefield types (spoiler,coupler,...)
//
#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchUtilities.h"

#include "AcceleratorModel/StdComponent/SWRFStructure.h"
#include "AcceleratorModel/StdComponent/TWRFStructure.h"
#include "AcceleratorModel/StdComponent/SectorBend.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/utils.h"
#include "IO/MerlinIO.h"
#include "TLAS/TLASimp.h"

#include <stdexcept>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>

#ifdef ENABLE_MPI
#include "MPI/merlin_mpi.hpp"

//Time for load leveling
//#include <ctime>
#include <time.h>
#endif

namespace {

#define COUT(x) cout<<std::setw(12)<<scientific<<std::setprecision(4)<<(x);

using namespace std;
using namespace PhysicalConstants;
using namespace PhysicalUnits;
using namespace ParticleTracking;

// needed to resolve gcc 3.2 ambiguity problem
//inline double pow(int x, int y) { return pow(double(x),double(y)); }

Point2D GetSliceCentroid(ParticleBunch::const_iterator first,
                         ParticleBunch::const_iterator last)
{
    Point2D c(0,0);
    double n=0;
    while(first!=last) {
        c.x += first->x();
        c.y += first->y();
        first++;
        n++;
    }
    return n>1 ? c/n : c;
}

PSvector GetSliceCentroid6D(ParticleBunch::const_iterator first,
                            ParticleBunch::const_iterator last)
{
    PSvector c(0);
    double n=0;
    while(first!=last) {
        c += *first;
        first++;
        n++;
    }
    if(n>1)
        c/=n;

    return c;
}

};



namespace ParticleTracking {

WakeFieldProcess::WakeFieldProcess (int prio, size_t nb, double ns, string aID)
        : ParticleBunchProcess(aID,prio),imploc(atExit),nbins(nb),nsig(ns),currentWake(0),Qd(),Qdp(),filter(0),
        wake_x(0),wake_y(0),wake_z(0),recalc(true),inc_tw(true),oldBunchLen(0)
{
    SetFilter(14,2,1);

    #ifdef ENABLE_MPI
    //Start the load leveling clock
//    t_initial = clock();
    clock_gettime(CLOCK_REALTIME, &t_initial);
    //int rank = MPI::COMM_WORLD.Get_rank();
    //cout << rank << "\t" << t_initial << endl << endl;
    #endif
}

WakeFieldProcess::~WakeFieldProcess()
{
    if(filter)
        delete filter;
}

size_t WakeFieldProcess::CalculateQdist()
{

    pair<double,double> v = currentBunch->GetMoments(ps_CT);
    double z0 = v.first;
    double sigz = v.second;

    // calculate binning ranges
//    cout << nsig << "\t" << sigz << "\t" << z0 << endl;
    zmin = -nsig*sigz+z0;
    zmax =  nsig*sigz+z0;
    dz = (zmax-zmin)/nbins;

    bunchSlices.clear();
    Qd.clear();
    Qdp.clear();

    // Qdp contains the slope of the charge distribution, smoothed using a filter
//    cout << "size 1: " << currentBunch->size() << endl;
    size_t lost = ParticleBinList(*currentBunch,zmin,zmax,nbins,bunchSlices,Qd,Qdp,filter);
//    cout << "size 2: " << currentBunch->size() << endl;
#ifndef NDEBUG
    ofstream os("qdist.dat");
    os<<zmin<<' '<<zmax<<' '<<dz<<endl;
    copy(Qd.begin(),Qd.end(),ostream_iterator<double>(os,"\n"));
#endif
//cout << "Rank: " << MPI::COMM_WORLD.Get_rank() << "\tQDIST: " << lost << endl;
    return lost;
}

// Smoothing filter takes the form of a set of coefficients
// calculated using the Savitzky-Golay technique
// n gives the width of the window on either side of the reference point
// m gives the order of the polynomial fitted to the points within the window
// d gives the order of the derivative required
// For CSR wake we need the first derivative
void WakeFieldProcess::SetFilter(int n, int m, int d)
{
	if(filter)
	{
		delete filter;
	}
	filter = new vector<double>;
	savgol(*filter, n, n, d, m);

	#ifndef NDEBUG
	ofstream os("filter.dat");
	copy(filter->begin(),filter->end(),ostream_iterator<double>(os,"\n"));
	#endif
}

void WakeFieldProcess::SetCurrentComponent (AcceleratorComponent& component)
{
    //	TWRFStructure* cavity = dynamic_cast<TWRFStructure*>(&component);
    //	WakePotentials* wake = cavity!=0 ? cavity->GetWakePotentials() : 0;

    WakePotentials* wake = component.GetWakePotentials();
    // if not initialize(=0) we assume that
    // WakeFieldProcess is responsible - for backward compatibility
    // in general expected process must be equal to this process
    if( wake &&
        wake->GetExpectedProcess()!=0 &&
        typeid(*(wake->GetExpectedProcess()))!=typeid(*this))
        wake=0;

    //if(wake!=0) cout<<GetID()<<endl;

	if(currentBunch!=0 && wake!=0)
	{
		clen = component.GetLength();
		switch(imploc)
		{
			case atCentre:
				impulse_s = clen/2.0;
				break;

			case atExit:
				impulse_s = clen;
				break;
		}

	        current_s = 0;
	        active = true;

		if(recalc || wake!=currentWake)
		{
			currentWake = wake;
			#ifndef ENABLE_MPI
			//If we are using an MPI build, we want to call Init after the master has been sent all the particles.
			Init();
			#endif
        	}


	}

	else
	{
		active = false;

		// check if we have a sector bend. If so, then
		// the bunch length will change and we need to rebin
		if(dynamic_cast<SectorBend*>(&component))
		{
			recalc = true;
		}
	}
}

void WakeFieldProcess::DoProcess(double ds)
{
	//Normal Code
	#ifndef ENABLE_MPI
//	cout << "WAKEFIELD PROCESS: " << currentBunch->size() << endl;
	current_s+=ds;
	if(fequal(current_s,impulse_s))
	{
		ApplyWakefield(clen);
		active = false;
	}
	#endif

	//MPI Code
	#ifdef ENABLE_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	//Step through the component
	current_s+=ds;

	//We only apply the wakefield if we are at the desired position: ImpulseLocation {atCentre,atExit};
	if(fequal(current_s,impulse_s))
	{
//	cout << "RANK: " << rank << "\tWAKEFIELD PROCESS::DO_PROCESS" << endl;

	//Define the MPI_Particle type as a structure
	//Should be easy to add extra fields to this, for macrocharge etc.
	const int coords = 6;
	int count = 1;
	int blockcounts[1];
	MPI::Aint offsets[1];
	MPI::Datatype oldtypes[1];
	offsets[0] = 0;
	oldtypes[0] = MPI_DOUBLE;
	blockcounts[0] = coords;

	//We now create the new datatype
	MPI::Datatype MPI_Particle = MPI::DOUBLE.Create_struct(count,blockcounts,offsets,oldtypes);
	//And commit it.
	MPI_Particle.Commit();

	MPI::Status merlin_mpi_status;

	//Before we hit the barrier, get how long it has taken to hit this point

	clock_gettime(CLOCK_REALTIME, &t_final);
	//timespec.tv_sec;
	//timespec.tv_nsec;
	//First find the number of seconds, convert to nanoseconds, then add the nanosecond difference.
	t_delta = ((t_final.tv_sec - t_initial.tv_sec)*1.0e9) + (t_final.tv_nsec - t_initial.tv_nsec);

	double time_per_particle = t_delta/currentBunch->size();	//In nanoseconds
	//cout << "time per particle: " << time_per_particle << endl;

	//We must collect the partial bunches from other nodes, then tell them to wait till the master has calculated the wakefield
	//Sync required first
	MPI::COMM_WORLD.Barrier();

	if(rank == 0)
	{
//		cout << "RANK 0 in WAKE" << endl;
		//First we must collect particles from all other "worker" nodes, and merge them in with the bunch on the "master" node.
		PSvector Particle_buffer;
		for(int n = 1; n<size; n++)
		{
                	//We probe the incomming send to see how many particles we are recieving.
	                MPI::COMM_WORLD.Probe(n, MPI_ANY_TAG, merlin_mpi_status);
	                //We find out from the probe the number of particles we have been sent
	                int recv_count = merlin_mpi_status.Get_count(MPI_Particle);
//	                cout << "Rank:" << rank << "\t Recv Count: " << recv_count << "\t From node: " << n << endl;

	                //We make a suitable buffer
	                double* particle_recv_buffer = new double[recv_count*coords];

        	        //We now do the Recv for real, 
	                MPI::COMM_WORLD.Recv(particle_recv_buffer, recv_count, MPI_Particle, n, 1, merlin_mpi_status);

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
	                        currentBunch->push_back(Particle_buffer);

/*
				//Debugging check
				
	                        for (int z = 0; z<6; z++)
				{
		                        if(isnan(Particle_buffer[z]))
		                        {
			                        //This is very bad so lets exit
			                        cout << "NAN recv from rank " << n << " in coord " << z+1 << " - particle number " << i << endl;
			                        MPI::COMM_WORLD.Abort(1);
	        	                }
	        	                else if(abs(Particle_buffer[z]) >= 1e5)
	        	                {
	        	                cout << "BIG" << endl;
	        	                cout << "NAN recv from rank " << n << " in coord " << z+1 << " - particle number " << i << endl;
	        	                MPI::COMM_WORLD.Abort(1);
	        	                }
				}//End check
*/				
        	        }

		delete [] particle_recv_buffer;

		}//Merge bunches now finished.
		
//		cout << endl <<"Master Bunch size: " << particle_count << endl << endl;
//		cout << "Pre-wake Charge: " << currentBunch->GetTotalCharge() << endl;
		//Now we have all the particles, Init() can be called.
		if(recalc)
		{
			Init();
		}

		//The "currentBunch" should now have all particles.
		//Now we can calculate the wakefield.

		ApplyWakefield(clen);
		active = false;

		//Now save the bunch size
		int particle_count = currentBunch->size();

		//Now we distribute the new bunch relative to the "load" on each worker node.
		double* load_balance = new double[size];

		//Fill the master node time
		load_balance[0] = time_per_particle;
		//cout << load_balance[0] << endl;
		//Recv the time values from other nodes
		for(int n = 1; n<size; n++)
		{
	                MPI::COMM_WORLD.Recv(&load_balance[n], 1, MPI::DOUBLE, n, 1, merlin_mpi_status);
//               		cout << load_balance[n] << endl;
		}

		//loop over the loads, then normalize.
		double total_load=0.0;
		for(int n = 0; n<size; n++)
		{
			total_load+=load_balance[n];
		}

		//Normalize
		double total_load_1=0.0;
		for(int n = 0; n<size; n++)
		{
			load_balance[n] = (total_load/load_balance[n]);
			total_load_1+=load_balance[n];
		}

		int total_particles_loaded = 0;
		//Fraction of total particles
		for(int n = 0; n<size; n++)
		{
			load_balance[n] = (int)(load_balance[n]*particle_count/total_load_1);
//			cout << "n: " << n << "\t npart: " << load_balance[n] << endl;
			total_particles_loaded+=(int)load_balance[n];
		}
//		cout << "Total particles:\t" << total_particles_loaded << endl;
//		cout << endl <<"Post wake size: " << particle_count << endl << endl;
//		cout << "Post-wake Charge: " << currentBunch->GetTotalCharge() << endl;

		//The bunch must now be convered into a format suitable for sending again.
		//The leftover extra particles can be added to the "local" machine bunch.
		//int particles_per_node = particle_count/size;
		//int remaining_particles = particle_count % size;
		int remaining_particles = particle_count - total_particles_loaded;

                double* particle_send_buffer = new double[particle_count*coords];
                for (int n=0; n<particle_count; n++)
                {
                        particle_send_buffer[(n*coords)+0] = currentBunch->GetParticles()[n].x();
                        particle_send_buffer[(n*coords)+1] = currentBunch->GetParticles()[n].xp();
                        particle_send_buffer[(n*coords)+2] = currentBunch->GetParticles()[n].y();
                        particle_send_buffer[(n*coords)+3] = currentBunch->GetParticles()[n].yp();
                        particle_send_buffer[(n*coords)+4] = currentBunch->GetParticles()[n].ct();
                        particle_send_buffer[(n*coords)+5] = currentBunch->GetParticles()[n].dp();
                }
		//We now have the full bunch in a particle buffer.
		//This must be sliced, and particles sent to each node

		//Since the macrocharge could have changed, we must also update the nodes with this new value.
		double macrocharge = currentBunch->GetTotalCharge()/particle_count;
		//cout << "MASTER MACROCHARGE: " << macrocharge << endl;

		//Resend
		//The old way:
		/*
		for(int n = 0; n<(size-1); n++)
                {
			//MPI::COMM_WORLD.Send(&particle_send_buffer[(particles_per_node*n*coord)+(particles_per_node+remaining_particles)], particles_per_node, MPI_Particle, (n+1), 1);
			MPI::COMM_WORLD.Send(&particle_send_buffer[(particles_per_node*n*coords)], particles_per_node, MPI_Particle, (n+1), 1);
		}
		*/

		//New adaptive resend
		int offset=0;
		for(int n = 1; n<size; n++)
                {
			MPI::COMM_WORLD.Send(&particle_send_buffer[offset], (int)load_balance[n], MPI_Particle, n, 1);
			offset+=(int)load_balance[n]*coords;
		}



		//Finally the local bunch must be reconstructed from the remaining particles.
		//First it must be cleared
		currentBunch->clear();

		//And then refilled.
                //for (int i = (particle_count-particles_per_node-remaining_particles); i<particle_count; i++)
                for (int i = (particle_count-(int)load_balance[0]-remaining_particles); i<particle_count; i++)
		//for (int i = 0; i<particles_per_node+remaining_particles; i++)
                {
                        Particle_buffer[0] = particle_send_buffer[(i*coords)+0];
                        Particle_buffer[1] = particle_send_buffer[(i*coords)+1];
                        Particle_buffer[2] = particle_send_buffer[(i*coords)+2];
                        Particle_buffer[3] = particle_send_buffer[(i*coords)+3];
                        Particle_buffer[4] = particle_send_buffer[(i*coords)+4];
                        Particle_buffer[5] = particle_send_buffer[(i*coords)+5];
                        //Push back each particle onto the cleared currentBunch
                        currentBunch->push_back(Particle_buffer);
                        
                        /*
                        //Debugging check
			for (int z = 0; z<6; z++)
			{
				if(isnan(Particle_buffer[z]))
				{
					//This is very bad so lets exit
					ofstream* badness = new ofstream("badness");
					currentBunch->Output(*badness);
					cout << "NAN particle at master " << z+1 << " - particle number " << i << endl;
					MPI::COMM_WORLD.Abort(1);
				}
			}//End check
			*/
                }
		delete [] particle_send_buffer;
		delete [] load_balance;

		//Now lets update the macrocharge on remote nodes.
		for(int n = 1; n<size; n++)
                {
			MPI::COMM_WORLD.Send(&macrocharge, 1, MPI::DOUBLE, n, 1);
		}

		//Sync
		MPI::COMM_WORLD.Barrier();
	}//End of rank 0 work.

	else
	{
		//Here we are a "worker" node, and will send particles to the master for the collective wakefield calculation
		//send bunches
		int particle_count = currentBunch->size();
	
		//We make an array to put the particle bunch data into.
		//It needs to be continuous in memory to work nicely with MPI
		//This is crude, but works.
		double* particle_send_buffer = new double[particle_count*coords];
		for (int n=0; n<particle_count; n++)
		{
			particle_send_buffer[(n*coords)+0] = currentBunch->GetParticles()[n].x();
			particle_send_buffer[(n*coords)+1] = currentBunch->GetParticles()[n].xp();
			particle_send_buffer[(n*coords)+2] = currentBunch->GetParticles()[n].y();
			particle_send_buffer[(n*coords)+3] = currentBunch->GetParticles()[n].yp();
			particle_send_buffer[(n*coords)+4] = currentBunch->GetParticles()[n].ct();
			particle_send_buffer[(n*coords)+5] = currentBunch->GetParticles()[n].dp();
		}
		
		//Send everything to the master node
//		cout << "RANK " << rank << ": SENDING" << endl;
		MPI::COMM_WORLD.Send(&particle_send_buffer[0], particle_count, MPI_Particle, 0, 1);

		//We now wait for wakefield to be calculated, and then recv some particles back.
		//Whilst the master node is busy calculating the Wakefield, we may as well use this time to clear the currentBunch.
		//This will be refilled with new particles.
		currentBunch->clear();

		//We then send the time per particle that it took to process the last batch to load level over the cluster.
		MPI::COMM_WORLD.Send(&time_per_particle, 1, MPI::DOUBLE, 0, 1);

		//We probe the incomming send to see how many particles we are recieving.
                MPI::COMM_WORLD.Probe(0, MPI_ANY_TAG, merlin_mpi_status);

		//We find out from the probe the number of particles we have been sent
                int recv_count = merlin_mpi_status.Get_count(MPI_Particle);
//		cout << "Rank:" << rank << "\t Recv Count: " << recv_count << endl;

		//We make a suitable buffer
		double* particle_recv_buffer = new double[recv_count*coords];

		//We now do the Recv for real, 
		MPI::COMM_WORLD.Recv(particle_recv_buffer, recv_count, MPI_Particle, 0, 1, merlin_mpi_status);

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
			//Particles[i][5] = particle_recv_buffer[(i*coords)+5];

			//Push back each particle onto the cleared currentBunch
			currentBunch->push_back(Particle_buffer);
			
			/*
			//Debugging check
			for (int z = 0; z<6; z++)
			{
				if(isnan(Particle_buffer[z]))
				{
					//This is very bad so lets exit
					cout << "NAN recv from master in coord " << z+1 << " - particle number " << i << endl;
					MPI::COMM_WORLD.Abort(1);
				}
			}//End check
			*/
		}
		
		//Finally we update the macrocharge per particle.
		double macrocharge;
		MPI::COMM_WORLD.Recv(&macrocharge, 1, MPI::DOUBLE, 0, 1, merlin_mpi_status);
		currentBunch->SetMacroParticleCharge(macrocharge);
		//cout << "Macrocharge: " << macrocharge << endl;

		current_s+=ds;
		delete [] particle_send_buffer;
		delete [] particle_recv_buffer;

		//Sync before carrying on with tracking
		MPI::COMM_WORLD.Barrier();
	}

	//And finally all processes must clear the MPI_Particle type.
	MPI_Particle.Free();

	//And now reset the intial clock
	//t_initial = clock();
	clock_gettime(CLOCK_REALTIME, &t_initial);

	}
//	cout << "RANK: " << rank << "\tWAKEFIELD PROCESS::END_PROCESS" << endl;
	#endif
}

void WakeFieldProcess::ApplyWakefield(double ds)
{
    // check if we are responsible for the current wakefield
    // a class derived from this class must 
    // include the appropriate check for its own wake potential type!
    if(typeid(WakePotentials*)!=typeid(currentWake)) return;

    // here we apply the wake field for
    // the step ds
    size_t n=0;
    double p0 = currentBunch->GetReferenceMomentum();

    // If the bunch length or binning has been changed,
    // we must recalculate the wakes
	// dk explicit check on bunch length
	if(recalc||oldBunchLen!=currentBunch->size())
        Init();

    // We always need to recalculate the transverse wake
    if(inc_tw)
        CalculateWakeT();

    // Now iterate over the bunch slices,
    // calculating the wakefield kicks using
    // linear interpolation between the values
    // at the slice boundaries

    double bload=0;

#define WAKE_GRADIENT(wake) (wake).empty() ? 0 : ((wake[nslice+1]-wake[nslice])/dz);

    double z=zmin;

    for(size_t nslice = 0; nslice<nbins; nslice++) {

        double gz = WAKE_GRADIENT(wake_z);
        double gx = WAKE_GRADIENT(wake_x);
        double gy = WAKE_GRADIENT(wake_y);


 

       for(ParticleBunch::iterator p=bunchSlices[nslice]; p!=bunchSlices[nslice+1]; p++) {
            double zz = p->ct()-z;
            double ddp = -ds*(wake_z[nslice]+gz*zz)/p0;
            p->dp() += ddp;
            bload += ddp;

            double dxp =  inc_tw? ds*(wake_x[nslice]+gx*zz)/p0 : 0;
            double dyp =  inc_tw? ds*(wake_y[nslice]+gy*zz)/p0 : 0;
           

            p->xp() = (p->xp()+dxp)/(1+ddp);
            p->yp() = (p->yp()+dyp)/(1+ddp);
        }
        z+=dz;
    }
    if(!currentWake->Is_CSR())
        currentBunch->AdjustRefMomentum(bload/currentBunch->size());
}

double WakeFieldProcess::GetMaxAllowedStepSize () const
{
    return impulse_s-current_s;
}

void WakeFieldProcess::Init()
{
    double Qt  = currentBunch->GetTotalCharge();

	//keep track of bunch length to be aware of modifications
	oldBunchLen=currentBunch->size();
    //cout << "INIT: CALCULATE QDIST: Total Charge: " << Qt << endl;
    size_t nloss = CalculateQdist();
    if(nloss!=0)
    {
        // Even though we have truncated particles, we still keep the
        // the bunch charge constant
        currentBunch->SetMacroParticleCharge(Qt/(currentBunch->size()));
        MerlinIO::warning()<<GetID()<<" (WakefieldProcess): "<<nloss<<" particles truncated"<<endl;
    }

    // Calculate the long. bunch wake.
    CalculateWakeL();
    recalc=false;

}



void WakeFieldProcess::CalculateWakeL()
{
    wake_z = vector<double>(bunchSlices.size(),0.0);
    double a0 = dz*fabs(currentBunch->GetTotalCharge())*ElectronCharge*Volt;

    // Estimate the bunch wake at the slice boundaries by
    // convolving the point-like wake over the current bunch
    // charge distribution Qd.
    //
    // Note that the distribution Qd is estimated at the
    // centre of each slice, not the slice boundary.
    //
    // CSR wake differs from classical wakes in two respects:
    // 1) wake from the tail of the bunch affects the head
    // 2) amplitude of wake depends on slope of charge distribution
    //    rather than directly on the distribution
    // Code to handle CSR wake added by A.Wolski 12/2/2003

    if(currentWake->Is_CSR())
    {
        for(int i=0; i<bunchSlices.size(); i++) {
            for(int j=1; j<i; j++) {
                wake_z[i] += Qdp[j]*(currentWake->Wlong((j-i+0.5)*dz))/dz;
            }
            wake_z[i]*=a0;
        }
    }
    else
    {
        for(int i=0; i<bunchSlices.size(); i++) {
            for(int j=i; j<bunchSlices.size()-1; j++) {
                wake_z[i] += Qd[j]*(currentWake->Wlong((j-i+0.5)*dz));
            }
            wake_z[i]*=a0;
        }
    }

#ifndef NDEBUG
    ofstream os("bunchWake.dat");
    os<<zmin<<'\t'<<zmax<<'\t'<<dz<<endl;
    copy(wake_z.begin(),wake_z.end(),ostream_iterator<double>(os,"\n"));
#endif

}

void WakeFieldProcess::CalculateWakeT()
{
    // This routine is rather cpu intensive.
    // First, calculate the transverse centroid of
    // each bunch slice by taking the mean of the
    // particle positions

    vector<Point2D> xyc;
    xyc.reserve(nbins);
    size_t i;
    for(i=0; i<nbins; i++)
        xyc.push_back(GetSliceCentroid(bunchSlices[i],bunchSlices[i+1]));

    // Now estimate the transverse bunch wake at the slice
    // boundaries in the same way we did for the longitudinal wake.

    double a0 = dz*(fabs(currentBunch->GetTotalCharge()))*ElectronCharge*Volt;
    wake_x = vector<double>(bunchSlices.size(),0.0);
    wake_y = vector<double>(bunchSlices.size(),0.0);
    for(i=0; i<bunchSlices.size(); i++)
    {
        for(int j=i; j<bunchSlices.size()-1; j++)
        {
            double wxy = Qd[j]*(currentWake->Wtrans((j-i+0.5)*dz));
            wake_x[i] += wxy*xyc[j].x;
            wake_y[i] += wxy*xyc[j].y;
        }
        wake_x[i]*=a0;
        wake_y[i]*=a0;
    }
}

void WakeFieldProcess::DumpSliceCentroids(ostream& os) const
{
    for(size_t i=0; i<nbins; i++) {
        os<<std::setw(4)<<i;
        os<<GetSliceCentroid6D(bunchSlices[i],bunchSlices[i+1]);
    }
}

void WakeFieldProcess::InitialiseProcess (Bunch& bunch)
{
    ParticleBunchProcess::InitialiseProcess(bunch);
    currentWake = 0;
    recalc = true;
}

// Calculate the Savitsky-Golay smoothing filter
// Adapted from Numerical Recipes in C
// This routine need only be executed once,
// when the WakeFieldProcess is initialized.
int powi(int i, int j)
{
	int p = 1;
	for(int m=0; m<j; m++)
		p *= i;
	return p;
};

void savgol(vector<double>& c, int nl, int nr, int ld, int m)
{
    Matrix<double> a(m+1,m+1);

    for(int i=0; i<=m; i++)
        for(int j=0; j<=i; j++)
        {
            double sum = 0.0;

            for(int k=-nl; k<=nr; k++)
                sum += powi(k, i) * powi(k, j);

            a(i,j) = sum;
            a(j,i) = sum;
        }

    vector<int> indx;
    double d;
    ludcmp(a,indx,d);

    Vector<double> b(m+1);
    for(int j=0; j<=m; j++)
        b(j) = 0.0;
    b(ld) = 1.0;

    lubksb(a,indx,b);

    c.clear();

    for(int k=-nl; k<=nr; k++)
    {
        double sum = b(0);
        double fac = 1.0;

        for(int mm=1; mm<=m; mm++)
            sum += b(mm) * (fac *= k);

        c.push_back(sum);
    }
}

}; // end namespace ParticleTracking

