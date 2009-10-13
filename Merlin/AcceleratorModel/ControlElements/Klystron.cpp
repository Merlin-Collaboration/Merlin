/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 13:55:32 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "AcceleratorModel/ControlElements/Klystron.h"
#include "AcceleratorModel/StdComponent/RFStructure.h"

using namespace std;

namespace {

	struct accumulate_vector_sum {
		accumulate_vector_sum() :z(0) {}
		void operator()(const RFStructure* cav) {
			z+= cav->GetVoltagePhasor();
		}
		Complex z;
	};

	struct set_cavity {
		set_cavity(double v, double phi) : v0(v), phi0(phi) {}

		void operator()(RFStructure* cav) {
			cav->SetVoltage(v0);
			cav->SetPhase(phi0);
		}

		const double v0;
		const double phi0;
	};

	struct scale_cavity {
		scale_cavity(const Complex& cs) : zs(cs) {}
		void operator()(RFStructure* cav) {
			Complex z = cav->GetVoltagePhasor();
			cav->SetVoltagePhasor(zs*z);
		}
		const Complex zs;
	};

}; // end anonumous namespace



Klystron::Klystron (const string& aName, const vector<RFStructure*>& cavs, Mode m)
	: ModelElement(aName),rf_cavs(cavs),kmode(m)
{
	Complex z0 = for_each(rf_cavs.begin(),rf_cavs.end(),accumulate_vector_sum()).z;
	V_k = abs(z0);
	phi_k =arg(z0);
	if(kmode==balanced) {
		const double v = V_k/rf_cavs.size();
		for_each(rf_cavs.begin(),rf_cavs.end(),set_cavity(v,phi_k));
	}
}

void Klystron::SetVoltage(double v)
{
	V_k = v;
	AdjustCavities();
}

void Klystron::SetPhase(double phi)
{
	phi_k = phi;
	AdjustCavities();
}

void Klystron::SetVoltagePhasor(const Complex& z)
{
	V_k = abs(z);
	phi_k = arg(z);
	AdjustCavities();
}

void Klystron::AdjustCavities()
{
	switch(kmode) {
		case balanced:
			for_each(rf_cavs.begin(),rf_cavs.end(),set_cavity(V_k/rf_cavs.size(),phi_k));
			break;
		case vector_sum:
			AdjustVectorSum(V_k,phi_k);
			break;
	}
}


const string& Klystron::GetType () const
{
	static const string typestr("Klystron");
	return typestr;
}

Klystron* Klystron::Copy () const
{
	return new Klystron(*this);
}

void Klystron::AppendBeamlineIndecies(std::vector<size_t>& ivec) const
{
	for(vector<RFStructure*>::const_iterator c = rf_cavs.begin(); c!=rf_cavs.end(); c++)
		ivec.push_back((*c)->GetBeamlineIndex());
}

void Klystron::AdjustVectorSum(double v, double phi)
{
	Complex vscale = Complex(v*cos(phi),v*sin(phi))/Complex(V_k*cos(phi_k),V_k*sin(phi_k));
	for_each(rf_cavs.begin(),rf_cavs.end(),scale_cavity(vscale));
	V_k = v;
	phi_k = phi;
}
