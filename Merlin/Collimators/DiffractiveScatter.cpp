/**
* Include for the Elastic scattering class
*/
#include "Collimators/DiffractiveScatter.h"

/**
* Include for the math headers - required for exp, sin, other mathmatical functions
*/
#include <cmath>

/**
* Include for io - std::cout etc
*/
#include <iostream>

/**
* Include for file output
*/
#include <fstream>

/**
* Include for the max() algorithm
*/
#include <algorithm>

/**
* Include for make_pair
*/
#include <utility>

/**
* Include for assorted numerial constants
*/
#include "NumericalUtils/NumericalConstants.h"

/**
* Include for assorted Physical constants
*/
#include "NumericalUtils/PhysicalConstants.h"

/**
* Include for assorted Physical units
*/
#include "NumericalUtils/PhysicalUnits.h"

/**
* Include for the random number generator
*/
#include "Random/RandomNG.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;

namespace ParticleTracking
{

/**
* Sets the minimum t value for generation
*/
void ppDiffractiveScatter::SetTMin(double tmin)
{
	t_min = tmin;
}

/**
* Sets the maximum t value for generation
*/
void ppDiffractiveScatter::SetTMax(double tmax)
{
	t_max = tmax;
}

/**
* Sets the step size in t for the differential cross section generation
*/
void ppDiffractiveScatter::SetTStepSize(double StepSize)
{
	t_step = StepSize;
}

/**
* Sets the minimum t value for generation
*/
void ppDiffractiveScatter::SetXiMin(double min)
{
	xi_min = min;
}

/**
* Sets the maximum t value for generation
*/
void ppDiffractiveScatter::SetXiMax(double max)
{
	xi_max = max;
}

/**
* Sets the step size in xi for the differential cross section generation
*/
void ppDiffractiveScatter::SetXiStepSize(double StepSize)
{
	xi_step = StepSize;
}

/**
* Gets the currently set minimum t value
*/
double ppDiffractiveScatter::GetTMin() const
{
	return t_min;
}

/**
* Gets the currently set maximum t value
*/
double ppDiffractiveScatter::GetTMax() const
{
	return t_max;
}

/**
* Gets the currently set t step size
*/
double ppDiffractiveScatter::GetTStepSize() const
{
	return t_step;
}

/**
* Gets the currently set minimum t value
*/
double ppDiffractiveScatter::GetXiMin() const
{
	return xi_min;
}

/**
* Gets the currently set maximum t value
*/
double ppDiffractiveScatter::GetXiMax() const
{
	return xi_max;
}

/**
* Gets the currently set xi step size
*/
double ppDiffractiveScatter::GetXiStepSize() const
{
	return xi_step;
}

/**
* Gets the Integrated elastic cross section
*/
double ppDiffractiveScatter::GetDiffractiveCrossSection() const
{
	return SigDiffractive;
}

/**
* Debug toggle - set to true or false to enable/disable debugging output
*/
void ppDiffractiveScatter::EnableDebug(bool flag)
{
	Debug = flag;
}

/**
* Generates the requried differential cross sections and integrates for the specified energy
*/
void ppDiffractiveScatter::GenerateDistribution(double energy)
{
	if(!Configured)
	{
		/*
		std::cout << "*******************************************************************************" << std::endl;
		std::cout << "*******   Generating pp Diffractive differential cross section   **************" << std::endl;
		std::cout << "*******************************************************************************" << std::endl;
		*/
		GenerateDsigDtDxi(energy);
		/*
		std::cout << "*******************************************************************************" << std::endl;
		std::cout << "*******************************************************************************" << std::endl;
		std::cout << "*************   Integrating differential cross section   **********************" << std::endl;
		std::cout << "*******************************************************************************" << std::endl;
		std::cout << "*******************************************************************************" << std::endl;
		IntegrateDsigDtDxi();

		std::cout << "*******************************************************************************" << std::endl;
		std::cout << "************   Diffractive Configuration generation done!   *******************" << std::endl;
		std::cout << "*******************************************************************************" << std::endl;
		*/
		Configured = true;
	}
}

ppDiffractiveScatter::~ppDiffractiveScatter()
{
} 
 

/**
* Generates the elastic differential cross section
* Places the results into the vectors t and DSig
* @param energy sqrt s
*/
/*
void ppDiffractiveScatter::GenerateDsigDtDxi(double energy)
{
	unsigned int nTSteps = (t_max - t_min) / t_step;
	unsigned int nXiSteps = (xi_max - xi_min) / xi_step;
	UniformT.reserve(nTSteps);
	UniformXi.reserve(nTSteps);
	DSig.reserve(nTSteps*nXiSteps);
	IntSig.reserve(nTSteps*nXiSteps);

	double s = (2 * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * energy) + (2 * pow(PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV,2));

	double sqrts = sqrt(s);
	std::cout << "Using " << nTSteps << " t bins and " << nXiSteps << " xi bins and sqrt s: " << sqrts << std::endl;

	if(!Debug)
	{
		for(unsigned int n = 0; n < nTSteps; n++)
		{
			UniformT.push_back((static_cast<double>(n) * nTSteps) + t_min);
			for(unsigned int m = 0; m < nXiSteps; m++)
			{
				UniformXi.push_back((static_cast<double>(m) * nXiSteps) + xi_min);
				DSig.push_back(PomeronScatter(UniformT[n],UniformXi[m],s));
			}
		}
	}
	else
	{
	}
}*/
void  __attribute__((optimize("O3,unsafe-math-optimizations"))) ppDiffractiveScatter::GenerateDsigDtDxi(const double energy)
{
	std::cout << "Call generateDsigDtDxi " << std::endl;
	const double s = (2 * PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV * energy) + (2 * pow(PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV,2));
	ss = s;
	std::cout << "s =" << s << std::endl;
	const int NN=10*N;
	//std::cout << "s =" << s << std::endl;
	t_step=(t_max-t_min)/NN;

	xi_step=(xi_max-xi_min)/NN;
	std::cout << "t_max" << "\t" << t_max << "\t" << "t_min" << "\t" << t_min << std::endl;
	std::cout << "xi_max" << "\t" << xi_max << "\t" << "xi_min" << "\t" << xi_min << std::endl;
	std::cout << "xi_step" << "\t" << xi_step << "\t" << "t_step" << "\t" << t_step << std::endl;
	double xdist[NN] = {0};
	double tdist[NN] = {0};
	//double xxx[NN];
	//double ttt[NN];
	//	double a[NN][NN];


	double sum=0;
//	std::cout << "PomeronScatter(-0.0001,0.0002,s) " << PomeronScatter(0.0001,0.0002,s) << std::endl;
//	std::cout << "PomeronScatter(-2,0.0003,s) " << PomeronScatter(2,0.0003,s) << std::endl;
//	std::cout << "PomeronScatter(-0.0001,0.002,s) " << PomeronScatter(0.0001,0.002,s) << std::endl;
//	std::cout << "PomeronScatter(-2,0.005,s) " << PomeronScatter(2,0.005,s) << std::endl;
	
	
	//std::cout << "dumping file" << std::endl;
	//std::ofstream *outfd = new std::ofstream("SD_Lookup_Table.dat");
	//(*outfd) << "#i" << "\t" << "xi" << "\t" << "t" <<  std::endl;
	
		
	
	
	//typedef std::pair <double, double> DoublePair;
	//std::vector <DoublePair> pairs;
	//double mrec=sqrt(s*x);
	for(int i=0;i<NN;i++)
	{
		const double t=t_min + i * t_step;
		//std::cout << "t " << t << std::endl;
		for(int j=0;j<NN;j++)
		{
			const double x=xi_min + j*xi_step;
			const double ds=PomeronScatter(t,x,s);
			sum=sum+ds;
			tdist[i]+=ds;
			xdist[j]+=ds;
			//std::cout << "x =" << x << "\tt =" << t << "\ttdist[" << i <<"] = "<< tdist[i] << "\txdist["<<j<<"] ="<< xdist[j] << "\tds =" << ds << std::endl;
			//unsigned int count = 0;
			//std::cout << "tdist\t" << tdist[i] << "xdist\t" << xdist[j] << std::endl;
			//std::cout << pairs[i] << "\t" <<  std::endl;
			//std::cout << t << "\t" << x << "\t" <<  ds << "\ttdist\t" << tdist[i] << "\txdist\t" << xdist[j]  << std::endl;
			//std::cout << "sum =" << sum << std::endl;
			 //if(j == NN-1)std::cout << "t" << t << "\tx" << x << "\tds " << ds << "\tsum " << sum << std::endl;
		}
		//if(i == NN-1) std::cout << "t" << t << "\tsum =" << sum << std::endl;
	}
//	for(int i=0;i<NN;i++)
//	{
//		std::cout << "tdist\t" << tdist[i] << "\txdist\t" << xdist[i] << std::endl;
//		std::cout << pairs[i].first << "\t" << pairs[i].second << std::endl;
//	}
	
	
	
	//outfd->close();
	//delete outfd;
	//std::cout << pairs[0].first << "\t" << pairs[0].second << std::endl;
	
	//for(int i=0;i<NN;i++)
//	{
//	std::cout << "tdist\t" << tdist[i] << "\txdist\t" << xdist[i] << std::endl;
//	}
	//sigma= 0.001*sum*xstep*tstep*1.6*pow(prop->A,0.33); // convert mbarn to barn
	//std::cout << "xi_step =" << xi_step << "\t" << "t_step" << t_step << std::endl;
	SigDiffractive= 0.001*sum * xi_step * t_step; // convert mbarn to barn  
	
	// convert histograms to normalised cumulants
	// Running total
	for(int i=1;i<NN;i++)
	{
		xdist[i]+=xdist[i-1];
		tdist[i]+=tdist[i-1];
	}
	
/*
	std::cout << "cumulant" << std::endl;
	for(int i=0;i<NN;i++)
	{
		std::cout << "tdist\t" << tdist[i] << "\txdist\t" << xdist[i] << std::endl;
	}
*/

	//Normalized
	for(int i=0;i<NN;i++)
	{
		xdist[i]/=xdist[NN-1];
		tdist[i]/=tdist[NN-1];
	}

/*
	std::cout << "normalized" << std::endl;
	for(int i=0;i<NN;i++)
	{
	std::cout << "tdist\t" << tdist[i] << "\txdist\t" << xdist[i] << std::endl;
	}
*/
	//Dump file
/*
	std::cout << "dumping file" << std::endl;
	std::ofstream *outfd = new std::ofstream("DiffractiveFile");
	std::ofstream *outfd1 = new std::ofstream("DiffractiveFile_lookup");
	(*outfd) << "#i" << "\t" << "xi" << "\t" << "t" <<  std::endl;
	for(int i=0;i<NN;i++)
	{
		(*outfd) << i << "\t" << xdist[i] << "\t" << tdist[i] << std::endl;;
	}
	outfd->close();
	delete outfd;
*/
	// convert to lookup tables
	int iseekt=0, iseekx=0;
	xarray[0]=tarray[0]=0;
	//int foo = std::cout.precision(10);
	//std::cout.precision(foo);
		
////	std::cout << "i\ttdist\ttdist+1\ttarray\tiseekt\tleftover\ttarget\tNN" << std::endl;
	for(int i=1;i<N;i++)
	{
		double target=double(i)/N;

		while(xdist[iseekx+1] < target)
		//while(xdist[iseekx] < target)
		{
			//std::cout << "xdist : " << xdist[iseekx] << "\t" << target << std::endl;
			iseekx++;
		}
		double leftover = target - xdist[iseekx];
		xarray[i] = (iseekx + std::min(1.0E0,leftover/(xdist[iseekx+1]-xdist[iseekx])))/NN;
		//std::cout << "Leftover: " << leftover << "\txarray:" << xarray[i] << "\tiseekx:" << iseekx << "\tMax:" << std::min(1.0E0,leftover/(xdist[iseekx+1]-xdist[iseekx])) << std::endl;
		while(tdist[iseekt+1] < target)
		{
			//std::cout << "tdist : " << tdist[iseekt] << "\t" << target << std::endl;
			iseekt++; 
		}
////		std::cout << iseekt << "\t" << tdist[iseekt+1] << "\t" << target << "\t" <<leftover/(tdist[iseekt+1]-tdist[iseekt]) << std::endl;
		leftover = target - tdist[iseekt];
		//tarray[i] = (iseekt + std::max(1.0E0,leftover/(tdist[iseekt+1]-tdist[iseekt])))/NN;
		tarray[i] = (static_cast<double>(iseekt) + std::min(1.0E0,leftover/(tdist[iseekt+1]-tdist[iseekt])))/static_cast<double>(NN);
		//std::cout << "Leftover: " << leftover << "\ttarray:" << tarray[i] << "\tiseekt:" << iseekt << "\tMax:" << std::min(1.0E0,leftover/(tdist[iseekt+1]-tdist[iseekt])) << std::endl;
		//tarray[1]=-tarray[1];
		//std::cout << "Leftover: " << leftover << "\ttarray:" << tarray[i] << std::endl;
//		std::cout << i << "\t" << tdist[iseekt] << "\t" << tdist[iseekt+1] << "\t" << tarray[i] << "\t" << iseekt << "\t" << leftover << "\t" <<
//		target << "\t" << static_cast<double>(iseekt) / static_cast<double>(NN) << std::endl;
	}





	//std::cout << "tarraysize\t" << (sizeof(tarray)/sizeof(*tarray))  << "\txarraysize\t" << (sizeof(xarray)/sizeof(*xarray))  << std::endl;

//	std::cout << "lookup table:" << std::endl;
//	for(int i=0;i<N;i++)
//	{
//		std::cout << "Generate i :  " << i << "\ttarray\t" << tarray[i] << "\txarray\t" << xarray[i] << "\tDtarray \t" << tarray[i+1]-tarray[i] << "\tDxarray \t" << xarray[i+1]-xarray[i] << std::endl;
//	}

	//outfd1->close();
	//delete outfd1;
	std::cout << "Nucleon Diffractive total cross section total "  << SigDiffractive * 1000.0 <<" mb" << std::endl;
	std::cout << "Sixtrack Diffractive total cross section total " << 0.00068*log(0.15*s) * 1000.0 <<" mb" << std::endl;
//	std::cout << "sum " << sum << " SigDiffractive " << SigDiffractive << std::endl;
	

}//End

/**
* Generates the SD differential cross section
* Places the results into the vectors t and DSig
*/
/*
void ppDiffractiveScatter::IntegrateDsigDtDxi()
{
	unsigned int nTSteps = UniformT.size();
	unsigned int nXiSteps = UniformXi.size();
	Sig.reserve(nTSteps*nXiSteps);
//	unsigned int seekt = 0;
	IntSig[0] = 0.0;
	std::cout << "INTEGRATE\t" << nTSteps << "\t" << nXiSteps << std::endl;
//	std::cout << Uniformt.size() << "\t" << DSig.size() << "\t" << IntSig.size() << std::endl;

	for(unsigned int n = 1; n < nTSteps*nXiSteps; n++)
	{
		double CurrentStepIntegral = (DSig[n] * nTSteps);//fix this
		SigDiffractive += CurrentStepIntegral;
		IntSig.push_back(SigDiffractive);
	}

	std::cout << "Diffractive Cross section (with peak): " << SigDiffractive * 1000 << std::endl;
	//Switch to normalized values to make our life easier
	for(unsigned int n = 1; n < nTSteps; n++)
	{
		for(unsigned int m = 1; m < nXiSteps; m++)
		{
			IntSig[n][m] /= SigDiffractive;
			//(*ofile) << UniformT[n] << "\t" << IntSig[n] << std::endl;
		}
	}

	InversionInterpolation = new Interpolation(IntSig, UniformT);
	
	for(unsigned int n=1; n <nSteps; n++)
	{
		double target = (static_cast<double>(n) / nSteps);
		
		try
		{
			sig_gen = (*InversionInterpolation)(target);
		}
		catch(Interpolation::BadRange& error)
		{
			std::cerr << "Bad Range in interpolation - requested: " << error.value << " but valid range is from " << error.valid_range.lower << " to "  << error.valid_range.upper << std::endl;
			std::cerr << "error in entry: " << n << " with total " << nSteps << std::endl;
			throw;
		}
		Sig.push_back(sig_gen);
	}

	LinearInterpolation = new Interpolation(Sig, 0, (1.0/nSteps));    // Interpolation of equally spaced data points
}
*/

/**
* Picks a t value from the generated distribution (including interpolation)
*/
double ppDiffractiveScatter::SelectT()
{

	double SigValue = RandomNG::uniform (0, 1.0);
	double t = (*LinearInterpolation)(SigValue);
	return t;
}

/**
* Picks an Xi value from the generated distribution (including interpolation)
*/
double ppDiffractiveScatter::SelectXi()
{
	double SigValue = RandomNG::uniform (0, 1.0);
	double xi = (*LinearInterpolation)(SigValue);
	return xi;
}

std::pair<double,double> ppDiffractiveScatter::Select()
{

	//static double fudge=1.0;
	
	//unsigned int count = 0;
	double xx,tt;
	static bool kilroy=false;
//	static ofstream ff;
	if(kilroy)
	{
		std::cout << "open file" << std::endl;
		kilroy=false;
//		ff.open("stuff.txt");
	}
	//for(int i=1;i<N;i++)
	//{
	//	std::cout << "Select   i:  " << i << "\ttarray\t" << tarray[i] << "\txarray\t" << xarray[i] << "\tDtarray \t" << tarray[i+1]-tarray[i] << "\tDxarray \t" << xarray[i+1]-xarray[i] << std::endl;
	//}
	retry:
	//count++;

	double rt = N*RandomNG::uniform(0,1);
	int it=int(rt); 
	double extra=rt-it;
	double deltat=0.0;
	if(it<(N-1))
	{
		deltat=tarray[it+1]-tarray[it];
		//std::cout << "\tdeltat : " << deltat << std::endl;
	}
	else deltat=1-tarray[it];
	tt=tarray[it] + extra * deltat;
	tt=t_min + tt *(t_max-t_min);
	//std::cout << "t = " << tt << "\t deltat = " << deltat << std::endl;
	if(tt < 0 || deltat < 0)
	{
		//std::cout << "Problem : it\t" << it << "  tarray["<< it+1 <<"]\t" << tarray[it+1] << "  tarray["<< it <<"]\t" << tarray[it] <<  "  deltat\t" << deltat << "  tt\t" << tt << std::endl;
		goto retry;
	}

	double rx = N*RandomNG::uniform(0,1);
	int ix=int(rx); 
	extra=rx-ix;
	double deltax=0.0;
	if(ix<(N-1))
	{
		deltax=xarray[ix+1]-xarray[ix];
	}
	else deltax=1-xarray[ix];
	xx=xarray[ix]+extra*deltax;
	xx=xi_min+xx*(xi_max-xi_min);
	//std::cout << "xx = " << xx << "\t deltax = " << deltax << std::endl;
	
	
	
	//double ds=PomeronScatter(tt,xx,s);
	//double sum = 77740400;
	//double sum = 101460;
	//int foo = std::cout.precision(10);
	//std::cout.precision(foo);
	double ds=PomeronScatter(tt,xx,ss)*0.001;
	//if(RandomNG::uniform(0,1)*sum > ds)
	//{
		//std::cout << "ds = " << ds  << std::endl;
	//	goto retry;
	//}
	//std::cout << "ds = " << ds <<std::endl;
	double ds2=SigDiffractive/(N*N*deltax*(xi_max-xi_min)*deltat*(t_max-t_min));
	//if(ds2<0)
	// std::cout << "Problem ds2 negative " << ds2 << "\tdeltax : " << deltax << "\tdeltat : " << deltat << std::endl;
	
	//
	//std::cout << "tt: " << tt << "/t xi : " << xx << "\tds: " << ds << "\tds2: " << ds2   << std::endl;
	//double P=ds/ds2;
	//if(P<RandomNG::uniform(0,1)) goto retry;
	
	static double fudge=1;
	double rat=fudge*ds/ds2;
	//double rat=fudge*ds2/ds;
	//std::cout << "rat :" << rat << std::endl;
	if(rat>1) fudge/=rat;//rat=fudge/rat;
	if(RandomNG::uniform(0,1)>rat)
	{
	//	std::cout <<"increase cut ratio by " << rat << "\tFudge: " << fudge << "\tds: " << ds << "\tds2: " << ds2 << std::endl;
	//	std::cout <<"SigDiffractive: " << SigDiffractive << "\tdeltax: " << deltax << "\tdeltat: " << deltat << "\tN: " << N << std::endl;
	//	std::cout << "rat: " << rat << "\t - " << ds << "\t - " << ds2 << "\t - " << fudge << std::endl;
	//	int foo = std::cout.precision(16);
	//	std::cout<< tarray[it-2] << "\t" << tarray[it-1] << "\t"  << tarray[it] << "\t" << tarray[it+1] << "\t" << "\t" << tarray[it+2] << it << std::endl;
	//	std::cout.precision(foo);
		//rat = fudge/rat;
		//fudge/=rat; //this is the original 
		//std::cout << "rejected : t = " << tt << "\tx = " << xx << "\tfudge = " <<  fudge << "\trat  " << rat << "\tds  " <<ds << "\tds2  " << ds2 << std::endl;
		goto retry;
	}
	//if (rat<RandomNG::uniform(0,1)) goto retry;
	//std::cout << "cycle : rat: " << rat << "\t  " << ds << "\t - " << ds2 << "\t - " << fudge << std::endl;
	//if (rat<RandomNG::uniform(0,1))
	//{
	 //std::cout << "t  " << tt << "\tx  " << xx << "\tds  " << ds << "\tds2  " << ds2 << std::endl;
	 //goto retry;
	//}
	//std::cout << "accepted : t = " << tt << "\tx = " << xx << "\t" << "ds  " <<ds << "\tds2: " << ds2 << std::endl;
	
	//if (rat * RandomNG::uniform(0,1)> 1) goto retry;

	//std::cout << "rat: " << rat << "\t - " << ds << "\t - " << ds2 << "\t - " << fudge << std::endl;
	


	//double tt = std::make_pair( t, sqrt(s*x)).first;

	//std::cout << count << std::endl;
	double mrec=sqrt(ss*xx);
	//std::cout << xx << "\t" << tt << std::endl;
//	ff << xx<<" "<<tt<<endl;
//	std::cout << "crash 3" << std::endl;
	return std::make_pair(tt,mrec);
	//return std::make_pair(tt,xx);
}

/**
* This is the high mass regge diffraction
*/

double ppDiffractiveScatter::PomeronScatter2(double tt, double x, double s)
{
	double sigma_xi = 0.1;
	double sigma_t = 0.1;
	//double t= - tt;
	double gausx = exp(-pow(x,2)/(2*sigma_xi));
	double gaust = exp(-pow(tt,2)/(2*sigma_t));
	//std::cout << gausx << "\t " << gaust << std::endl;
	return pow(x,2);
	
}
// 
inline double  __attribute__((optimize("O3,unsafe-math-optimizations"))) __attribute__ ((hot)) ppDiffractiveScatter::PomeronScatter(const double tt, const double x, const double s)   const
{
	//std::cout << "Call PomeronScatter SD " << std::endl;
	const double t= - tt;

	// Parameters for the resonance term in the background
	const double Mproton = 0.938272013;
	const double Mpion = 0.1349766;
	const double Mmin2 = pow(Mproton+Mpion,2);
	const double ml01 = 1.44;  //mass of the resonance P11 D13 G15 F17
        const double ml02 = 1.52;
        const double ml03 = 1.68;
	const double ml04 = 2.19;
	const double GammaL1 = 0.325; // width of the resonance
	const double GammaL2 = 0.13;
	const double GammaL3 = 0.14;
	const double GammaL4 = 0.45;
	const double cl01 = 3.07;   // coupling coefficient from the data fit
        const double cl02 = 0.4149; 
	const double cl03 = 1.108  ;
	const double cl04 = 0.9515;
	const double Mcut = 3; // this is chosen from the fit on the crosss section data
	const double xi_c = pow(Mcut,2)/s;
	const double xi_th = Mmin2/s; // (M_p + M_pion)^2/s
	const double Mmin2bar = pow(Mproton-Mpion,2);
	if(x <= xi_th){return 0;}

	double cc[4][3];
/*	cc[0][0] = 0.881148;	//A
	cc[0][1] = 3.94056;	//B
	cc[0][2] = 0.0220505;		//C

	//ppr
	cc[1][0] = 2.42997;
	cc[1][1] = 3.11514;
	cc[1][2] = 0.104746;
	
	//rrp
	cc[2][0] = 6.28648;
	cc[2][1] = 4.05376;
	cc[2][2] = 8.63609;

	//rrr
	cc[3][0] = 167.618;
	cc[3][1] = 11.5978;
	cc[3][2] = 54.256849;
*/
	cc[0][0] = 0.624529;	//A
	cc[0][1] = 2.5835;	//B
	cc[0][2] = 0;		//C

	//ppr
	cc[1][0] = 3.09088;
	cc[1][1] = 4.51487;
	cc[1][2] = 0.186211;
	
	//rrp
	cc[2][0] = 4.0;
	cc[2][1] = 3.03392;
	cc[2][2] = 10.0;

	//rrr
	cc[3][0] = 177.217;
	cc[3][1] = 5.86474;
	cc[3][2] = 21.0029;



	const double q = sqrt((x*s-Mmin2)*(x*s-Mmin2bar)/(4*x*s));
	const double ql1 = sqrt( (pow(ml01,2) - Mmin2) * (pow(ml01,2) - Mmin2bar) /(4*pow(ml01,2)) );
	const double ql2 = sqrt( (pow(ml02,2) - Mmin2) * (pow(ml02,2) - Mmin2bar) /(4*pow(ml02,2)) );
	const double ql3 = sqrt( (pow(ml03,2) - Mmin2) * (pow(ml03,2) - Mmin2bar) /(4*pow(ml03,2)) );
	const double ql4 = sqrt( (pow(ml04,2) - Mmin2) * (pow(ml04,2) - Mmin2bar) /(4*pow(ml04,2)) );
	//std::cout << ql1 << '\t' << ql2 << '\t' << ql3 <<  '\t' << ql4 << std::endl;
	const double gammaL01 = GammaL1*pow(q/ql1,3)*((1 + 5*ql1)/(1 + 5*q));
	const double gammaL02 = GammaL2*pow(q/ql2,5)*pow(((1 + 5*ql2)/(1 + 5*q)),2);
	const double gammaL03 = GammaL3*pow(q/ql3,7)*pow(((1 + 5*ql3)/(1 + 5*q)),3);
	const double gammaL04 = GammaL4*pow(q/ql4,9)*pow(((1 + 5*ql4)/(1 + 5*q)),4);

	const double R = ( (cl01/x) * (ml01*gammaL01) / ( pow( (x*s - pow(ml01,2) ),2) + pow(ml01*gammaL01,2)) 
		    +(cl02/x) * (ml02*gammaL02) / ( pow( (x*s - pow(ml02,2) ),2) + pow(ml02*gammaL02,2))
		    +(cl03/x) * (ml03*gammaL03) / ( pow( (x*s - pow(ml03,2) ),2) + pow(ml03*gammaL03,2))
		    +(cl04/x) * (ml04*gammaL04) / ( pow( (x*s - pow(ml04,2) ),2) + pow(ml04*gammaL04,2)) ) 
		    *exp(13.5*(t + 0.05));// Normalization factors Sandy's note 
		     //* sqrt(565/s)*exp(13.5*(t + 0.05));// Normalization factors Sandy's note 
//	double BRMatch = - 588.20982975 *exp(13.5*(t + 0.05))*(x - xi_th)/(xi_c - xi_th);
	const double BRMatch = -  ( (cl01/xi_c) * (ml01*gammaL01) / ( pow( (xi_c*s - pow(ml01,2) ),2) + pow(ml01*gammaL01,2))
		    +(cl02/xi_c) * (ml02*gammaL02) / ( pow( (xi_c*s - pow(ml02,2) ),2) + pow(ml02*gammaL02,2))
		    +(cl03/xi_c) * (ml03*gammaL03) / ( pow( (xi_c*s - pow(ml03,2) ),2) + pow(ml03*gammaL03,2))
		    +(cl04/xi_c) * (ml04*gammaL04) / ( pow( (xi_c*s - pow(ml04,2) ),2) + pow(ml04*gammaL04,2)) )
		    *exp(13.5*(t + 0.05)) *(x - xi_th)/(xi_c - xi_th);


//			   *sqrt(565/s)*exp(13.5*(t + 0.05))*(x - xi_th)/(xi_c - xi_th);

	//std::cout << "t = " << t << std::endl;	
	if(t > -0.25)
	{
		//std::cout << "t less than 1.15" << std::endl;
		if(x > xi_th && x <= xi_c)
		{
			const double Axi_c = (0.4+0.5*t)*pow(s,0.08) * pow(xi_c,-1.08 -0.5*t)	//ppp
			+(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(xi_c,-1.6125 -0.5*t)	//ppr
			+(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * pow(xi_c,-0.015 - 1.86*t)	//rrp
			+(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(xi_c,-0.5475 - 1.86*t)	//rrr
			+1.14591559 * pow(3.52142 - 2.79*t,2) * pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))	//pion
			* fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor
		
						
			const double Aprimexi_c =(0.4+0.5*t)*pow(s,0.08) * (-1.08 - 0.5*t) *pow(xi_c,-2.08 - 0.5*t)	//ppp
			+(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * (-1.6125 - 0.5*t) *pow(xi_c,-2.6125 -0.5*t)	//ppr
			+(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * (-0.015 - 1.86*t) *pow(xi_c,-1.015 - 1.86*t)	//rrp
			+(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * (-0.5475 - 1.86*t) * pow(xi_c,-1.5475 - 1.86*t)	//rrr
			+1.14591559 * pow(3.52142 - 2.79*t,2) * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2) 
			* ((1 - 1.86 * (-0.0182185 + t)) *pow(xi_c,-1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))
			+ (pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*(-0.4525)*pow(s*xi_c,-1.4525)+13.63*0.0808*pow(s*xi_c,0.0808-1))));	//form factor
	 	
			const double d = ((xi_c-xi_th)*Aprimexi_c-Axi_c)/pow(xi_c-xi_th,2);
			const double e = Aprimexi_c -2*((xi_c-xi_th)*Aprimexi_c-Axi_c)/(xi_c-xi_th);
			
			const double B = d * pow(x - xi_th,2) + e * (x - xi_th);
			//return B + R + BRMatch;
			//std::cout << t << "\t" << x << "\t" <<  B << "\t" << R << "\t" << B+R+BRMatch <<std::endl;
			return B + R + BRMatch; 
			//std::cout << t << "\t" << x << "\t" <<  BTot << std::endl;
			//return BTot;
			
		}
		else
		{
			return (0.4+0.5*t)*pow(s,0.08) * pow(x,-1.08 -0.5*t)	//ppp
			+(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(x,-1.6125 -0.5*t)	//ppr
			+(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) 	* pow(x,-0.015 - 1.86*t)	//rrp
			+(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(x,-0.5475 - 1.86*t)	//rrr
			+1.14591559 * pow(3.52142 - 2.79*t,2) * pow(x,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*x,-0.4525) + 13.63 *pow(s*x,0.0808))	//pion
			* fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor
			//std::cout << "t" << "\t" << "x" << "\t" <<  "A" << std::endl;
			//std::cout << t << "\t" << x << "\t" <<  A << std::endl;
			//return A;
		}
	
	}else if(t > -1.15)
	{
		//std::cout << "t less than 1.15" << std::endl;
		if(x > xi_th && x <= xi_c)
		{
			const double Axi_c = (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * pow(xi_c,-1.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			+(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(xi_c,-1.6125 -0.5*t)	//ppr
			+(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * pow(xi_c,-0.015 - 1.86*t)	//rrp
			+(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(xi_c,-0.5475 - 1.86*t)	//rrr
			+1.14591559 * pow(3.52142 - 2.79*t,2) * pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))	//pion
			* fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor
		
						
			const double Aprimexi_c = (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * (-1.08 - 0.5*t) *pow(xi_c,-2.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			+(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * (-1.6125 - 0.5*t) *pow(xi_c,-2.6125 -0.5*t)	//ppr
			+(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * (-0.015 - 1.86*t) *pow(xi_c,-1.015 - 1.86*t)	//rrp
			+(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * (-0.5475 - 1.86*t) * pow(xi_c,-1.5475 - 1.86*t)	//rrr
			+1.14591559 * pow(3.52142 - 2.79*t,2) * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2)
			* ((1 - 1.86 * (-0.0182185 + t)) *pow(xi_c,-1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))
			+ (pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*(-0.4525)*pow(s*xi_c,-1.4525)+13.63*0.0808*pow(s*xi_c,0.0808-1))));	//form factor
	 	
			const double d = ((xi_c-xi_th)*Aprimexi_c-Axi_c)/pow(xi_c-xi_th,2);
			const double e = Aprimexi_c -2*((xi_c-xi_th)*Aprimexi_c-Axi_c)/(xi_c-xi_th);
			
			const double B = d * pow(x - xi_th,2) + e * (x - xi_th);
			//return B + R + BRMatch;
			//std::cout << t << "\t" << x << "\t" <<  B << "\t" << R << "\t" << B+R+BRMatch <<std::endl;
			return B + R + BRMatch; 
			//std::cout << t << "\t" << x << "\t" <<  BTot << std::endl;
			//return BTot;
			
		}
		else
		{
			return (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * pow(x,-1.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			+(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(x,-1.6125 -0.5*t)	//ppr
			+(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) 	* pow(x,-0.015 - 1.86*t)	//rrp
			+(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(x,-0.5475 - 1.86*t)	//rrr
			+1.14591559 * pow(3.52142 - 2.79*t,2) * pow(x,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*x,-0.4525) + 13.63 *pow(s*x,0.0808))	//pion
			* fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor
			//std::cout << "t" << "\t" << "x" << "\t" <<  "A" << std::endl;
			//std::cout << t << "\t" << x << "\t" <<  A << std::endl;
			//return A;
		}
	
	}
	else
	{	//std::cout << "t bigger or equal than 1.15" << std::endl;
		if(x > xi_th && x <= xi_c)
		{
			const double Axi_c = (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * pow(xi_c,-1.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			*(1 + 0.4597*(fabs(t)-1.15) + 5.7575 * pow((fabs(t)-1.15),2))
			+(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(xi_c,-1.6125 -0.5*t)	//ppr
			+(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * pow(xi_c,-0.015 - 1.86*t)	//rrp
			+(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(xi_c,-0.5475 - 1.86*t)	//rrr
			+1.14591559 * pow(3.52142 - 2.79*t,2) * pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808))	//pion
			* fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor
		
						
			const double Aprimexi_c = (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * (-1.08 - 0.5*t) *pow(xi_c,-2.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			*(1 + 0.4597*(fabs(t)-1.15) + 5.7575 * pow((fabs(t)-1.15),2))
			+(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * (-1.6125 - 0.5*t) *pow(xi_c,-2.6125 -0.5*t)	//ppr
			+(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) * (-0.015 - 1.86*t) *pow(xi_c,-1.015 - 1.86*t)	//rrp
			+(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * (-0.5475 - 1.86*t) * pow(xi_c,-1.5475 - 1.86*t)	//rrr
			+1.14591559 * pow(3.52142 - 2.79*t,2) * fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2) 
			* ((1 - 1.86 * (-0.0182185 + t)) *pow(xi_c,-1.86 * (-0.0182185 + t)) * (31.79*pow(s*xi_c,-0.4525) + 13.63 *pow(s*xi_c,0.0808)) 
			+ (pow(xi_c,1 - 1.86 * (-0.0182185 + t)) * (31.79*(-0.4525)*pow(s*xi_c,-1.4525)+13.63*0.0808*pow(s*xi_c,0.0808-1))));	//form factor
	 	
			const double d = ((xi_c-xi_th)*Aprimexi_c-Axi_c)/pow(xi_c-xi_th,2);
			const double e = Aprimexi_c -2*((xi_c-xi_th)*Aprimexi_c-Axi_c)/(xi_c-xi_th);
			
			const double B = d * pow(x - xi_th,2) + e * (x - xi_th);
			//return B + R + BRMatch;
			//std::cout << t << "\t" << x << "\t" <<  B << "\t" << R << "\t" << B+R+BRMatch <<std::endl;
			return B + R + BRMatch; 
			//std::cout << t << "\t" << x << "\t" <<  BTot << std::endl;
			//return BTot;
			
		}
		else
		{
			return (cc[0][0]*exp(cc[0][1]*t) + cc[0][2]) * pow(s,0.08) * pow(x,-1.08 - 0.5*t) * (t/(t - 0.05))	//ppp
			*(1 + 0.4597*(fabs(t)-1.15) + 5.7575 * pow((fabs(t)-1.15),2))
			+(cc[1][0]*exp(cc[1][1]*t) + cc[1][2]) * pow(s,-0.4525) * pow(x,-1.6125 -0.5*t)	//ppr
			+(cc[2][0]*exp(cc[2][1]*t) + cc[2][2]) * pow(s,0.08) 	* pow(x,-0.015 - 1.86*t)	//rrp
			+(cc[3][0]*exp(cc[3][1]*t) + cc[3][2]) * pow(s,-0.4525) * pow(x,-0.5475 - 1.86*t)	//rrr
			+1.14591559 * pow(3.52142 - 2.79*t,2) * pow(x,1 - 1.86 * (-0.0182185 + t)) * (31.79*pow(s*x,-0.4525) + 13.63 *pow(s*x,0.0808))	//pion
			* fabs(t) * pow(1 - 1.40845*t,-4) * pow(3.52142 -t,-2) * pow(-0.0182185 + t,-2);	//form factor
			//std::cout << "t" << "\t" << "x" << "\t" <<  "A" << std::endl;
			//std::cout << t << "\t" << x << "\t" <<  A << std::endl;
			//return A;
		}
		
	}
}
}//End namespace ParticleTracking
