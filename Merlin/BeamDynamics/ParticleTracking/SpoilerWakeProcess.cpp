/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Created: June, 2006
//
/////////////////////////////////////////////////////////////////////////
//
// additional spoiler wakefield code thanks to Roger Barlow and Adriana Bunglau
// see EUROTeV Report 2006-051 for the physics 
// the present implementation is different
// Modified by D.Kruecker 18.2.2008
//

#include "merlin_config.h"
#include "AcceleratorModel/SpoilerWakePotentials.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchUtilities.h"
#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"
#include "BeamDynamics/ParticleTracking/SpoilerWakeProcess.h"
#include "AcceleratorModel/Aperture.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <vector>

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

namespace
{
  inline long double powd (double x, int y)
  { long double product (1);
    for (int i=1; i<=y; i++)
      product *= x;
    return product;
  }
}

namespace ParticleTracking{

// Constructor

SpoilerWakeProcess::SpoilerWakeProcess(int modes, int prio, size_t nb, double ns)
                   : WakeFieldProcess (prio, nb, ns,"SPOILERWAKE") {  
                   if(modes>5) {
                   	cout<<"SpoilerWakeProcess: Max number of modes is 5 - Set to 5!"<<endl;
                   	modes=5;
                   }
                   nmodes = modes;
}

// Destructor

SpoilerWakeProcess:: ~SpoilerWakeProcess(){}
  

// Calculates the moments Cm for each slice 

double SpoilerWakeProcess::CalculateCm(int mode, int slice)
{	
          double x=0;
          for(ParticleBunch::iterator p=bunchSlices[slice]; p!=bunchSlices[slice+1]; p++)
          {
          double r = sqrt(powd(p->x(),2)+powd(p->y(),2));
          double theta = atan2(p->y(),p->x());
          x += powd(r,mode)*cos(mode*theta);
          }
          return x;
}


// Calculates  the moments Sm for each slice 

double SpoilerWakeProcess::CalculateSm(int mode, int slice)
{
          double x=0;
          for(ParticleBunch::iterator p=bunchSlices[slice]; p!=bunchSlices[slice+1]; p++)
          {
          double r = sqrt(powd(p->x(),2)+powd(p->y(),2));
          double theta = atan2(p->y(),p->x());
          x += powd(r,mode)*sin(mode*theta);
          }
          return x;
}

// Calculate the transverse wake with modes 

void SpoilerWakeProcess::CalculateWakeT(double dz, int currmode)
 {     
            vector <double> w(nbins); 
    for (size_t slice=0; slice<nbins; slice++)
          {
             w[slice] = spoiler_wake-> Wtrans(slice*dz, currmode);
             
           } 
  
     for (size_t slice=0; slice<nbins; slice++)
     { 	 
            int i = slice;
          wake_ct[currmode][i]=0;
          wake_st[currmode][i]=0;

             for(size_t j=slice; j<bunchSlices.size()-1; j++)
              {
		wake_ct[currmode][i] +=  w[j-i]*Cm[currmode][j];
	        wake_st[currmode][i] += w[j-i]*Sm[currmode][j];  
              }  
           }	
}

// This function calculates the longitudinal wake
// with modes

void SpoilerWakeProcess::CalculateWakeL(double dz, int currmode)
{  
    vector<double> w(nbins);    
    for (size_t slice=0; slice<nbins; slice++)
          {
             w[slice] = spoiler_wake-> Wlong(slice*dz, currmode);
         
            }   
     for (size_t i=0; i<nbins; i++)
      { 
          wake_cl[currmode][i]=0;
          wake_sl[currmode][i]=0;
          for(size_t j=i; j<bunchSlices.size()-1; j++)
          {
	    wake_cl[currmode][i] += w[j-i]*Cm[currmode][j];
	    wake_sl[currmode][i] += w[j-i]*Sm[currmode][j];
          }  
        }	
}

void SpoilerWakeProcess::ApplyWakefield(double ds){

	// down casting to its own wake potential type
	// we are not responsible if it is a different type 
	spoiler_wake=dynamic_cast<SpoilerWakePotentials*>(currentWake);
	if(spoiler_wake==0) return;

	// If the bunch length or binning has been changed,
    // we must recalculate the wakes
    if(recalc||oldBunchLen!=currentBunch->size())
        Init();

	{int m;
	for(m=1;m<=nmodes;m++){     
	    for (size_t  n=0;n<nbins;n++){
	    	 Cm[m][n]=CalculateCm(m,n);
		 Sm[m][n]=CalculateSm(m,n);
	    }
        }}
	double wake_x,wake_y,wake_z;          
	double macrocharge=currentBunch->GetTotalCharge()/currentBunch->size();
        
	double a0 = macrocharge*ElectronCharge*Volt;
	a0 *= 3E8;

	double p0 = currentBunch->GetReferenceMomentum();
	//dk ? if(recalc) Init();
	double bload=0;
	       
#define WAKE_GRADIENT(wake) ((wake[currmode][nslice+1]-wake[currmode][nslice])/dz); 
  
	double z=zmin;
	
	for(int currmode=1; currmode<=nmodes; currmode++){

		CalculateWakeT(dz, currmode);	
		CalculateWakeL(dz, currmode);	

		for(size_t nslice = 0; nslice<nbins; nslice++) {
			double g_ct = WAKE_GRADIENT(wake_ct);
			double g_st = WAKE_GRADIENT(wake_st);
			double g_cl = WAKE_GRADIENT(wake_cl);
			double g_sl = WAKE_GRADIENT(wake_sl);
			for(ParticleBunch::iterator p=bunchSlices[nslice]; p!=bunchSlices[nslice+1]; p++)
			{
				double r = sqrt (powd(p->x(),2) + powd(p->y(),2)); 
				double theta = atan2(p->y(),p->x());
				double zz = p->ct()-z;
				double wxc = cos((currmode-1)*theta)*(wake_ct[currmode][nslice]+g_ct*zz); 
				double wxs = sin((currmode-1)*theta)*(wake_st[currmode][nslice]+g_st*zz);
				double wys = cos((currmode-1)*theta)*(wake_st[currmode][nslice]+g_st*zz);
				double wyc = sin((currmode-1)*theta)*(wake_ct[currmode][nslice]+g_ct*zz);
				wake_x = currmode*powd(r,currmode-1)*(wxc+wxs);
				wake_y = currmode*powd(r,currmode-1)*(wys-wyc);
				wake_x*=a0;
				wake_y*=a0;
				double wzc = cos(currmode*theta)*(wake_cl[currmode][nslice]+g_cl*zz);
				double wzs = sin(currmode*theta)*(wake_sl[currmode][nslice]+g_sl*zz);
				wake_z = powd(r,currmode)*(wzc-wzs);
				wake_z*= a0;
				double ddp = -ds*wake_z/p0;
				p->dp() += ddp;
				bload += ddp;
				double dxp = inc_tw? ds*wake_x/p0 : 0;
				double dyp = inc_tw? ds*wake_y/p0 : 0;
				p->xp() = (p->xp()+dxp)/(1+ddp);
				p->yp() = (p->yp()+dyp)/(1+ddp);
			}
			z+=dz;  
    		}
	}   
}
}; //end namespace ParticleTracking

