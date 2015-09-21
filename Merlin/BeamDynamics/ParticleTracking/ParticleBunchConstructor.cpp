/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/04/14 16:01:43 $
// $Revision: 1.9 $
//
/////////////////////////////////////////////////////////////////////////

#include "BasicTransport/NormalTransform.h"
#include <memory>
// RandomNG
#include "Random/RandomNG.h"
// ParticleBunchConstructor
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"
// Include for pi constant.
#include "NumericalUtils/NumericalConstants.h"

namespace ParticleTracking {

inline double RandomGauss(double variance, double cutoff)
{
    return cutoff==0 ? RandomNG::normal(0,variance) :  RandomNG::normal(0,variance,cutoff);
}

ParticleBunchConstructor::ParticleBunchConstructor (const BeamData& beam, size_t npart, DistributionType dist)
        : np(npart),dtype(dist),cutoffs(0),beamdat(beam),itsFilter(0),M(NormalTransform(beam)),force_c(false)
{}

ParticleBunchConstructor::~ParticleBunchConstructor ()
{
    if(itsFilter)
        delete itsFilter;
}

void ParticleBunchConstructor::SetBunchData (const BeamData& beam)
{
    beamdat = beam;
    M.R = NormalTransform(beam);
}

void ParticleBunchConstructor::SetNumParticles (size_t npart)
{
    assert(npart>0);
    np=npart;
}

void ParticleBunchConstructor::SetDistributionCutoff (double cut)
{
    cutoffs = PSvector(fabs(cut));
}

void ParticleBunchConstructor::SetDistributionCutoff (const PSvector& cut)
{
    cutoffs=cut;
}

//Bunch* ParticleBunchConstructor::ConstructBunchDistribution (int bunchIndex) const
void ParticleBunchConstructor::ConstructBunchDistribution (int bunchIndex) const
{
    PSvector p;

//    PSvectorArray pbunch = pbunch1;
    // First we generate npart particles in "normalised" phase
    // space, after which we transform them to "real" phase
    // space using M.

    // The first particle is *always* the centroid particle
    double dp2 = pow(beamdat.sig_dp,2);
    double dz2 = pow(beamdat.sig_z,2);
    double rx,ry;
    double u;

    p.x()=beamdat.x0;
    p.xp()=beamdat.xp0;
    p.y()=beamdat.y0;
    p.yp()=beamdat.yp0;
    p.dp()=0;
    p.ct()=beamdat.ct0;
    p.type() = -1.0;
    p.location() = -1.0;
    p.id() = 0;
    p.sd() = 0.0;
    pbunch.push_back(p);

    size_t i;

    PSvector xm = p; // used for calculating mean

    switch(dtype) {
    case normalDistribution:
        for(i=1; i<np;) {
            p.x()	= RandomGauss(beamdat.emit_x,cutoffs.x());
            p.xp()	= RandomGauss(beamdat.emit_x,cutoffs.xp());
            p.y()	= RandomGauss(beamdat.emit_y,cutoffs.y());
            p.yp()	= RandomGauss(beamdat.emit_y,cutoffs.yp());
            p.dp()	= RandomGauss(dp2,cutoffs.dp());
            p.ct()	= RandomGauss(dz2,cutoffs.ct());

            M.Apply(p);
            p+=pbunch.front(); // add centroid
            p.type() = -1.0;
            p.location() = -1.0;
            p.id() = 0;
            p.sd() = 0.0;

            if(itsFilter==0 || itsFilter->Apply(p)) {
                pbunch.push_back(p);
                xm += p;
                i++;
            }

        }
        if(force_c) {
            xm/=np;
            xm-=pbunch.front();
            PSvectorArray::iterator pp=pbunch.begin();
            pp++;
            for(;pp!=pbunch.end(); pp++)
                (*pp)-=xm;
        }
        break;
    case flatDistribution:
        rx = sqrt(beamdat.emit_x);
        ry = sqrt(beamdat.emit_y);
        for(i=1; i<np;) {
            p.x()	= RandomNG::uniform(-rx,rx);
            p.xp()	= RandomNG::uniform(-rx,rx);
            p.y()	= RandomNG::uniform(-ry,ry);
            p.yp()	= RandomNG::uniform(-ry,ry);
            p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
            p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
            M.Apply(p);
            p+=pbunch.front(); // add centroid
            p.type() = -1.0;
            p.location() = -1.0;
            p.id() = 0;
            p.sd() = 0.0;
            

            if(itsFilter==0 || itsFilter->Apply(p)) {
                pbunch.push_back(p);
                i++;
            }
        }
        break;
        case skewHaloDistribution:
        case ringDistribution:
	           rx = sqrt(beamdat.emit_x);
	           ry = sqrt(beamdat.emit_y);		  
//	           pbunch.pop_back();
	           for(i=1; i<np;) {
	               u = RandomNG::uniform(-pi,pi);
	               p.x()	= rx * cos(u);
	               p.xp()	= rx * sin(u);
	               u = RandomNG::uniform(-pi,pi);
	               p.y()	= ry * cos(u);
	               p.yp()	= ry * sin(u);
	               p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
	               p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
	               M.Apply(p);
	               p+=pbunch.front(); // add centroid
            p.type() = -1.0;
            p.location() = -1.0;
            p.id() = 0;
            p.sd() = 0.0;
			//cout<<p<<endl;
	               if(itsFilter==0 || itsFilter->Apply(p)) {
	                   pbunch.push_back(p);
	                   i++;
	               }
	           }
        break;
        case horizontalHaloDistribution1:
	           rx = sqrt(beamdat.emit_x);
	           ry = sqrt(beamdat.emit_y);
	           for(i=1; i<np;) {
	               u = RandomNG::uniform(-pi,pi);
	               p.x()	= rx * cos(u);
	               p.xp()	= rx * sin(u);
	               u = RandomNG::uniform(-pi,pi);
	               p.y()	= 0.0;
	               p.yp()	= 0.0;
	               p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
	               p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
                       // cout<<p<<endl;
		       M.Apply(p);
	               p+=pbunch.front(); // add centroid
	               p.type() = -1.0;
	               p.location() = -1.0;
	               p.sd() = 0.0;
	               p.id() = i;
	               if(itsFilter==0 || itsFilter->Apply(p)) {
	                   pbunch.push_back(p);
	                   i++;
	               }
	           }
	break;
        case verticalHaloDistribution1:
	           rx = sqrt(beamdat.emit_x);
	           ry = sqrt(beamdat.emit_y);
	           for(i=1; i<np;) {
	               u = RandomNG::uniform(-pi,pi);
	               p.x()	= 0.0;
	               p.xp()	= 0.0;
	               u = RandomNG::uniform(-pi,pi);
	               p.y()	= ry * cos(u);
	               p.yp()	= ry * sin(u);
	               p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
	               p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
	              cout<<"rx\t"<< rx << endl;
		      cout<<"beamdat.emit_x"<<'\t'<<beamdat.emit_x<<endl;
                       M.Apply(p);
	               p+=pbunch.front(); // add centroid
	               p.type() = -1.0;
	               p.location() = -1.0;
	               p.sd() = 0.0;
	               p.id() = i;
	               //cout<<p<<endl;
	               if(itsFilter==0 || itsFilter->Apply(p)) {
	                   pbunch.push_back(p);
	                   i++;
	               }
	           }
	break;
        case horizontalHaloDistribution2:
	           rx = sqrt(beamdat.emit_x);
	           ry = sqrt(beamdat.emit_y);
	           for(i=1; i<np;) {
	               u = RandomNG::uniform(-pi,pi);
	               p.x()	= rx * cos(u);
	               p.xp()	= rx * sin(u);
	               u = RandomNG::uniform(-pi,pi);
	               p.y()	= RandomGauss(beamdat.emit_y,cutoffs.y());
	               p.yp()	= RandomGauss(beamdat.emit_y,cutoffs.yp());
	               p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
	               p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);	               
		      //cout << p << endl;
		       M.Apply(p);
	               p+=pbunch.front(); // add centroid
	               p.type() = -1.0;
	               p.location() = -1.0;
	               p.sd() = 0.0;
	               p.id() = i;
	               if(itsFilter==0 || itsFilter->Apply(p)) {
	                   pbunch.push_back(p);
	                   i++;
	               }
	           }
	break;
	case horizontalHaloDistribution3:
	           rx = sqrt(beamdat.emit_x);
	           ry = sqrt(beamdat.emit_y);
	           for(i=1; i<np;) {
	               u = RandomNG::uniform(-pi,pi);
			if(i%2 == 0)
	               p.x()	= 0.00151881506107592;// half jaw width + 1 microm (TCP.D6L7.B1)
			if(i%2 != 0)
	               p.x()	= -0.00151881506107592;// half jaw width + 1 microm (TCP.D6L7.B1)
	               p.xp()	= 0;
	               u = RandomNG::uniform(-pi,pi);
	               p.y()	= 0;//RandomGauss(beamdat.emit_y,cutoffs.y());
	               p.yp()	= 0;//RandomGauss(beamdat.emit_y,cutoffs.yp());
	               p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
	               p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);	               
		      //cout << p << endl;
		       //M.Apply(p);
	               p+=pbunch.front(); // add centroid
	               p.type() = -1.0;
	               p.location() = -1.0;
	               p.sd() = 0.0;
	               p.id() = i;
	               if(itsFilter==0 || itsFilter->Apply(p)) {
	                   pbunch.push_back(p);
	                   i++;
	               }
	           }
	break;
        case verticalHaloDistribution2:
	           rx = sqrt(beamdat.emit_x);
	           ry = sqrt(beamdat.emit_y);
	           for(i=1; i<np;) {
	               u = RandomNG::uniform(-pi,pi);
	               p.x()	= RandomGauss(beamdat.emit_x,cutoffs.x());
	               p.xp()	= RandomGauss(beamdat.emit_x,cutoffs.xp());
	               u = RandomNG::uniform(-pi,pi);
	               p.y()	= ry * cos(u);
	               p.yp()	= ry * sin(u);
	               p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
	               p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
	               M.Apply(p);
	               p+=pbunch.front(); // add centroid
	               p.type() = -1.0;
	               p.location() = -1.0;
	               p.sd() = 0.0;
	               p.id() = i;
	               if(itsFilter==0 || itsFilter->Apply(p)) {
	                   pbunch.push_back(p);
	                   i++;
	               }
	           }
	break;
	case tuneTestDistribution:
		{
		// from m to n sigma
		rx = sqrt(beamdat.emit_x);
		ry = sqrt(beamdat.emit_y);
		//Nominal LHC		
		//~ double sigx = 267.067E-6; //TCP.C6L7		
		//~ double sigx =  0.0002966510271229613; //HEL
		
		//HEL (HL)
		//~ double sigx = 0.000258007;	//HEL
		double sigx = 0.000266382;	//TCP.C6L7
		
		// with M.apply (TCP)
		//~ double first = 1.1E-4; //5.2 sigma
		//~ double last = 1.3E-4; //6 sigma
		//~ double first = 0.7E-5; //~1 sigma
		//~ double last = 1.3E-4; //~6 sigma
		
		//without M.apply			
		double first = 7*sigx; 
		double last = 9*sigx; 	
		//~ double first = 1*sigx; 
		//~ double last = 10*sigx; 	
			
		double partsi = np-2;
		double steps = (last-first)/partsi;
		double intpart;
			
		for(i=1; i<np;) {
			u = RandomNG::uniform(-1,1);

			//~ double test = rx*u;	
			p.x()	= (first + ((i-1) * steps));
			p.xp() = 0;			
						
			//~ p.x()	= rx*u;			
			//~ double fracpart = modf((p.x()/rx), &intpart);
			//~ fracpart = sqrt(pow(fracpart,2));
			//~ p.xp() = rx *  sin( acos(fracpart) );
			
			//~ p.xp() = rx *  sin( acos(p.x()/sqrt(rx)) );
			//~ p.xp() = -2.4E-5;
			//~ p.xp() = -p.x()*0.016;

			//~ cout << "\n\t\t\tParticleBunchConstructor: x = " << p.x() << " xp = " << p.xp() << endl; 
		
			//~ p.xp()	= rx *  sin( sqrt((p.x()/rx)*(p.x()/rx)) *pi);
			//~ p.xp()	= rx *  (cos(  sqrt(pow(p.x(),2) / pow(rx, 2)) * 0.5 *pi ));
			//~ p.xp()	= sqrt( pow(p.x(),2) + pow(rx, 2));
			//~ p.xp() = rx-p.x();
			//~ u = RandomNG::uniform(-1,1);
			//~ if (u>=0){
				//~ p.xp() *= -1;
			//~ }
		
			p.y()	= 0.0;
			p.yp()	= 0.0;
			p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
			p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);	               
			//~ M.Apply(p);
			p+=pbunch.front(); // add centroid

			if(itsFilter==0 || itsFilter->Apply(p))  {
				pbunch.push_back(p);
				i++;
			}
		}
	}
	break;
    };

    //return new ParticleBunch(beamdat.p0,beamdat.charge,pbunch);
}

Bunch* ParticleBunchConstructor::ConstructBunch (int bunchIndex) const
{
	ParticleBunchConstructor::ConstructBunchDistribution();
	return new ParticleBunch(beamdat.p0,beamdat.charge,pbunch);
}

void ParticleBunchConstructor::ForceCentroid (bool fc)
{
    force_c = fc;
}

ParticleBunchFilter::~ParticleBunchFilter ()
{
    // Nothing to do
}

} //end namespace ParticleTracking
