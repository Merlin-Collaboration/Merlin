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

#include "merlin_config.h"
#include "Collimators/SpoilerWakePotentials.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchUtilities.h"
#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"
#include "Collimators/SpoilerWakeProcess.h"
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
	inline double powd (double x, int y)
	{
		double product (1);
		for (int i=1; i<=y; i++)
		{
			product *= x;
		}
		return product;
	}
}

namespace ParticleTracking{

// Constructor

SpoilerWakeProcess::SpoilerWakeProcess(int modes, int prio, size_t nb, double ns)
  : WakeFieldProcess (prio, nb, ns)
{
	// cout<<" SpoilerWakeProcess "<<modes<<endl;
	nmodes = modes;
	Cm=new double*[modes+1];
	Sm=new double*[modes+1];
	wake_ct=new double*[modes+1];
	wake_st=new double*[modes+1];
	wake_cl=new double*[modes+1];
	wake_sl=new double*[modes+1];
	for (int i=0;i<=nmodes;i++)
	{
		Cm[i]=new double[1000];
		Sm[i]=new double[1000];
		wake_ct[i]=new double[1000];
		wake_st[i]=new double[1000];
		wake_cl[i]=new double[1000];
		wake_sl[i]=new double[1000];
	}
}

// Destructor

SpoilerWakeProcess:: ~SpoilerWakeProcess()
{
	for (int i=0;i<nmodes;i++)
	{
		delete[] Cm[i];
		delete[] Sm[i];
		delete[] wake_ct[i];
		delete[] wake_st[i];
		delete[] wake_cl[i];
		delete[] wake_sl[i];
	}

	delete[]Cm;
	delete[]Sm;
	delete[] wake_ct;
	delete[] wake_st;
	delete[] wake_cl;
	delete[] wake_sl;
}
  

// Calculates the moments Cm for each slice 
double SpoilerWakeProcess::CalculateCm(int mode, int slice)
{	
	// cout<<"Cm  mode and slice "<<mode<<" "<<slice<<endl;
	double x=0;
	for(ParticleBunch::iterator p=bunchSlices[slice]; p!=bunchSlices[slice+1]; p++)
	{
		double r = sqrt(powd(p->x(),2)+powd(p->y(),2));
		double theta = atan2(p->y(),p->x());
		x += powd(r,mode)*cos(mode*theta);
		// cout<<" value x "<<x<<endl;
	}
	return x;
}


// Calculates  the moments Sm for each slice 
double SpoilerWakeProcess::CalculateSm(int mode, int slice)
{
	double x=0;
	// cout<<"Sm  mode and slice "<<mode<<" "<<slice<<endl;
	for(ParticleBunch::iterator p=bunchSlices[slice]; p!=bunchSlices[slice+1]; p++)
	{
		double r = sqrt(powd(p->x(),2)+powd(p->y(),2));
		double theta = atan2(p->y(),p->x());
		x += powd(r,mode)*sin(mode*theta);
		// cout<<" value x "<<x<<endl;
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
		//cout<<" Wake "<<currmode<<" "<<slice<<" "<<w[slice]<<endl;
		//cout << w[slice] << "\t" << slice << "\t" << dz << endl;
	}

	for (size_t slice=0; slice<nbins; slice++)
	{
		int i = slice;
		wake_ct[currmode][i]=0;
		wake_st[currmode][i]=0;

		for(size_t j=slice; j<bunchSlices.size()-1; j++)
		{
			wake_ct[currmode][i] += w[j-i]*Cm[currmode][j];
		        wake_st[currmode][i] += w[j-i]*Sm[currmode][j];
		}
	}
}

// This function calculates the longitudinal wake with modes
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

void SpoilerWakeProcess::ApplyWakefield(double ds)//  int nmodes)
{
	//CalculateQdist();

	//pair<double,double> v=currentBunch ->GetMoments(ps_CT);
	//double z0=v.first;
	//double sigz=v.second;
	//ofstream *logfile = new ofstream("Output/logfile.dat");
	//logfile->precision(12);
	//ofstream *bunchafterfile = new ofstream("Output/wakefield.dat");
	//currentBunch->Output(*bunchafterfile);
	//delete bunchafterfile;
	
	spoiler_wake=(SpoilerWakePotentials*) currentWake;
	for(int m=1;m<=nmodes;m++)
	{
		for (size_t  n=0;n<nbins;n++)
		{
			Cm[m][n]=CalculateCm(m,n);
			Sm[m][n]=CalculateSm(m,n);
		}
	}

	double wake_x,wake_y,wake_z;
	double macrocharge=currentBunch->GetTotalCharge()/currentBunch->size();
	double a0 = macrocharge*ElectronCharge*Volt;
	a0 /= 4*pi*FreeSpacePermittivity;
	double p0 = currentBunch->GetReferenceMomentum();

	if(recalc) Init();
	double bload=0;

	#define WAKE_GRADIENT(wake) ((wake[currmode][nslice+1]-wake[currmode][nslice])/dz);

	for(int currmode=1; currmode<=nmodes; currmode++)
	{
		CalculateWakeT(dz, currmode);
		CalculateWakeL(dz, currmode);
		double z=zmin;
		int iparticle=0;

		for(size_t nslice = 0; nslice<nbins; nslice++)
		{
			// cout<<" trig stuff "<<currmode<<" "<<nslice<<" "<<Cm[currmode][nslice]<<" "<<Sm[currmode][nslice]<<endl;
			double g_ct = WAKE_GRADIENT(wake_ct);
			double g_st = WAKE_GRADIENT(wake_st);
			double g_cl = WAKE_GRADIENT(wake_cl);
			double g_sl = WAKE_GRADIENT(wake_sl);
			g_ct=g_st=g_cl=g_sl=0;
			int number_particles=0;
			for(ParticleBunch::iterator p=bunchSlices[nslice]; p!=bunchSlices[nslice+1]; p++)
			{
				number_particles++;
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
		//	*logfile<<currmode<<"\t"<<r<< "\t" << theta << "\t" << wake_ct[currmode][nslice] <<"\t"<<wake_st[currmode][nslice]<<"\t" << wxc << "\t" << wxs << "\t" << wyc << "\t" << wys << endl;    //<<p->x()<<"\t"<<p->y()<<endl;
				double wzc = cos(currmode*theta)*(wake_cl[currmode][nslice]+g_cl*zz);
				double wzs = sin(currmode*theta)*(wake_sl[currmode][nslice]+g_sl*zz);
				wake_z = powd(r,currmode)*(wzc-wzs);
				wake_z*= a0;
				double ddp = -wake_z/p0;
				p->dp() += ddp;
				bload += ddp;
				double dxp = inc_tw? wake_x/p0 : 0;
				double dyp = inc_tw? wake_y/p0 : 0;
				p->xp() = (p->xp()+dxp)/(1+ddp);
				double oldpy=p->yp();
				p->yp() = (p->yp()+dyp)/(1+ddp);
			// *logfile <<" new... slice "<<nslice<<" particle "<<(++iparticle)<< "\tangles "<< p->xp()<<"\t"<< p->yp()<<"\t"<<dxp<<"\t"<<dyp<<endl;
			}
			//*logfile << "slice number: " << nslice << "\tslice size: " << number_particles << endl;			
			z+=dz;

		//	ostringstream bunchsave;
		//	bunchsave << "Output/bunch_" << nslice << ".bunch";
		//        ofstream *bunchsavefile = new ofstream(bunchsave.str().c_str());
		 //       currentBunch->Output(*bunchsavefile);
		//	bunchsavefile->close();
		//        delete bunchsavefile;
		}
	}
//	ofstream *bunchafterfile1 = new ofstream("Output/wakefield1.dat");
//	currentBunch->Output(*bunchafterfile1);
//	delete bunchafterfile1;
//	delete logfile;
}
}; //end namespace ParticleTracking
