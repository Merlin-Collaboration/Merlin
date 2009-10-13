#include "HaloTracker.h"
#include "merlin_config.h"
#include <iostream> 
#include <fstream>
#include <cmath>

#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
#include "BeamDynamics/ParticleTracking/CollimateParticleProcess.h"
#include "NumericalUtils/NumericalConstants.h"

#include "Random/RandomNG.h"
#include "QuadIntegrator.h"

using namespace ParticleTracking;

namespace {
    
    // local function to calculate a halo as used by STRUCT
    ParticleBunch* STRUCT_Bunch(const BeamData& bdat, double nx1,
        double nx2, double ny1, double ny2, double dpp,
        int npart);
}; // end anonymous namespace

HaloTracker::HaloTracker(AcceleratorModel::Beamline beamline, const BeamData& beamdat)
: collimate_halo(true),dump_particles(true),use_struct_dist(false),
n1x(0),n2x(0),n1y(0),n2y(0),x(0),xp(0),y(0),yp(0),dp(0),bdat(beamdat),bline(beamline)
{}

void HaloTracker::Run(size_t npart)
{
    // Construct flat halo bunch
    ParticleBunch* bunch0;
    
    if(!use_struct_dist) {
        bunch0 = new ParticleBunch(bdat.p0);
        while(npart--) {
            Particle p(0);
            p.x()  = RandomNG::uniform(-x,x)+bdat.x0;
            p.xp() = RandomNG::uniform(-xp,xp)+bdat.xp0;
            p.y()  = RandomNG::uniform(-y,y)+bdat.y0;
            p.yp() = RandomNG::uniform(-yp,yp)+bdat.yp0;
            p.dp() = RandomNG::uniform(-dp,dp);
            p.ct() = 0;
            bunch0->AddParticle(p);
        }
    }
    else
        bunch0 = STRUCT_Bunch(bdat,n1x,n2x,n1y,n2y,dp,npart);
    
    // Set up the tracker
    ParticleTracker tracker(bline,bunch0,true);
    
    // Here we override the library quadrupole
    // integrator with our own local version
    // (see QuadIntgrator.h)
    tracker.RegisterIntegrator(new QuadIntegrator);

    ofstream lossSummaryOS("data/loss_summary.dat");
    
    if(collimate_halo) {
        
        // Loss summary file is data/loss_summary.dat
        CollimateParticleProcess* collproc = 
            new CollimateParticleProcess(2,COLL_AT_EXIT,&lossSummaryOS);
        
        // Create loss particle files in data/
        collproc->CreateParticleLossFiles(true,"data/");
        
        // Abort tracking when 100% particles are lost.
        collproc->SetLossThreshold(100.0);
        
        // Index the loss particle files
        collproc->IndexParticles(true);

		// Turn on/off spoiler scattering
		collproc->ScatterAtSpoiler(scatter_at_spoiler);
        
        tracker.AddProcess(collproc);
    }
    
    // When 100% of the particles are lost, the process
    // throws a MerlinException. Thus we must wrap
    // tracker.Run() in a try block.
    
    try {
        tracker.Run();
        if(dump_particles) {
            ofstream os("data/particles.out.dat");
            tracker.GetTrackedBunch().Output(os);
        }
    }
    catch(MerlinException& err) {
        cout<<err.Msg()<<endl;
    }	
}

void HaloTracker::SetSTRUCTlimits(double nX1, double nX2, double nY1, 
                                  double nY2, double dpp)
{
    use_struct_dist=true;
    n1x=nX1;
    n2x=nX2;
    n1y=nY1;
    n2y=nY2;
    dp=dpp;
}


// set the halo limits (meter, radian)
void HaloTracker::SetHaloLimits(double x1, double xp1, double y1, double yp1, double dp1)
{
    x=x1;
    xp=xp1;
    y=y1;
    yp=yp1;
    dp=dp1;
}


// set the halo limits normalised to the nominal sigma
void HaloTracker::SetHaloLimitsN(double nx, double nxp, double ny, double nyp, double dp1)
{
    double gamma_x = (1+pow(bdat.alpha_x,2))/bdat.beta_x;
    double gamma_y = (1+pow(bdat.alpha_y,2))/bdat.beta_y;
    
    double sigx = sqrt(bdat.emit_x*bdat.beta_x);
    double sigxp = sqrt(bdat.emit_x*gamma_x);
    double sigy = sqrt(bdat.emit_y*bdat.beta_y);
    double sigyp = sqrt(bdat.emit_y*gamma_y);
    
    x = nx*sigx;
    xp = nxp*sigxp;
    y = ny*sigy;
    yp=nyp*sigyp;
    dp = dp1;
}


namespace {

    ParticleBunch* STRUCT_Bunch(const BeamData& bdat, double nx1,
        double nx2, double ny1, double ny2, double ddp, 
        int npart)
    {
        ParticleBunch* bunch = new ParticleBunch(bdat.p0);
        
        double gamma_x = (1+pow(bdat.alpha_x,2))/bdat.beta_x;
        double gamma_y = (1+pow(bdat.alpha_y,2))/bdat.beta_y;
        
        double sigx = sqrt(bdat.emit_x*bdat.beta_x);
        double sigxp = sqrt(bdat.emit_x*gamma_x);
        double sigy = sqrt(bdat.emit_y*bdat.beta_y);
        double sigyp = sqrt(bdat.emit_y*gamma_y);
        
        double logNx = log(nx2/nx1);
        double logNy = log(ny2/ny1);
        
        while(npart>0) {
            
            Particle p;
            
            // generate 1/r distribution
            double Ax = nx1*exp(logNx*RandomNG::uniform(0,1));
            double Ay = ny1*exp(logNy*RandomNG::uniform(0,1));
            
            double phix = RandomNG::uniform(0,twoPi);
            double phiy = RandomNG::uniform(0,twoPi);
            
            p.x() = sigx*Ax*cos(phix);
            p.xp()= sigxp*Ax*sin(phix);
            p.y() = sigy*Ay*cos(phiy);
            p.yp()= sigyp*Ay*sin(phiy);
            p.dp()= RandomNG::normal(0,ddp*ddp);
            
            bunch->AddParticle(p);
            npart--;
        }
        
        return bunch;
    }
    
}; 
