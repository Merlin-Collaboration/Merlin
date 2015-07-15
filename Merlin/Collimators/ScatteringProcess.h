#ifndef _h_ScatteringProcess
#define _h_ScatteringProcess 1

#include "merlin_config.h"

#include "BeamModel/PSvector.h"

#include "Collimators/Material.h"


/*

Definition of the virtual ScatteringProcess class and
also several child classes derived from it
 
Created RJB 23 October 2012
Modified HR 2013

*/

class ScatteringProcess {
public:
	double sigma; // integrated cross section for this process

protected:
	double E0;    // Beam energy
	Material* mat; // material MaterialProperties of the collimator being hit
public:
	// The first function must be provided for all child classes, and probably the second as well
	virtual bool scatter(PSvector& p, double E)=0;
	virtual void Configure(Material* matin, double E0in){mat=matin; E0=E0in;}
	virtual std::string GetProcessType() const {return "ScatteringProcess";}
};

class Rutherford:public ScatteringProcess
{
	double tmin;
public:
	void Configure(Material* matin, double E0in);
	bool scatter(PSvector& p, double E);
	std::string GetProcessType() const{return "Rutherford";}
};

class NucleonElasticSimple:public ScatteringProcess
{
	double slope; 
public:
	void Configure(Material* mat, double P0);
	bool scatter(PSvector& p, double E);
	std::string GetProcessType() const {return "NucleonElasticSimple";}
};

class NucleusElastic:public ScatteringProcess
{
	double slope; 
public:
	void Configure(Material* mat, double P0);
	bool scatter(PSvector& p, double E);
	std::string GetProcessType() const {return "NucleusElastic";}
};

class NucleonDiffractiveSimple:public ScatteringProcess
{
	double slope,Mx_lo2,Mx_hi2;
	double s; 
public:
	void Configure(Material* mat, double P0);
	bool scatter(PSvector& p, double E);
	std::string GetProcessType() const {return "NucleonDiffractiveSimple";}
};

class Inelastic:public ScatteringProcess
{
public:
	void Configure(Material* mat, double P0);
	bool scatter(PSvector& p, double E){return false;} // Particle inevitably lost
	std::string GetProcessType() const {return "Inelastic";}
};

//HR Sixtrack-like scattering

class SixNucleonElastic : public ScatteringProcess
{
    double slope;
public:
    void Configure(Material* mat, double P0);
	bool scatter(PSvector& p, double E);
	std::string GetProcessType() const {return "SixtrackNucleonElastic";}

};


class SixSingleDiffractive : public ScatteringProcess
{
    double b;
    double xm2;
    double comsqd;
public:
    void Configure(Material* mat, double P0);
	bool scatter(PSvector& p, double E);
	std::string GetProcessType() const {return "SixtrackSingleDiffractive";}

};

#endif
