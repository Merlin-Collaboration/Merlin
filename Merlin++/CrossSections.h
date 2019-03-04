/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_CrossSections
#define _h_CrossSections 1

#include <iostream>
#include <cmath>

#include "merlin_config.h"

#include "PSvector.h"

#include "Material.h"
#include "DiffractiveScatter.h"
#include "ElasticScatter.h"

#include "utils.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

namespace Collimation
{

class CrossSections
{
public:

	//default constructor
	CrossSections();

	//overloaded constructor
	CrossSections(Material* mat, double E, int scattertype);

	~CrossSections();

	inline bool operator==(const CrossSections& rhs)
	{
		if((this->sig_pN_tot_ref != rhs.sig_pN_tot_ref) || (this->scat_type != rhs.scat_type) || (this->lambda_tot !=
			rhs.lambda_tot))
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}

	void ConfigureCrossSections(double E0);

	double GetTotalMeanFreePath();

	void Set_E0(double a)
	{
		E0 = a;
	}
	double Get_E0() const
	{
		return E0;
	}

	void Set_sig_pN_tot_ref(double a)
	{
		sig_pN_tot_ref = a;
	}
	double Get_sig_pN_tot_ref() const
	{
		return sig_pN_tot_ref;
	}

	void Set_sig_pN_inel_ref(double a)
	{
		sig_pN_inel_ref = a;
	}
	double Get_sig_pN_inel_ref() const
	{
		return sig_pN_inel_ref;
	}

	void Set_sig_R_ref(double a)
	{
		sig_R_ref = a;
	}
	double Get_sig_R_ref() const
	{
		return sig_R_ref;
	}

	void Set_sig_R(double a)
	{
		sig_R = a;
	}
	double Get_sig_R() const
	{
		return sig_R;
	}

	void Set_sig_pp_tot(double a)
	{
		sig_pp_tot = a;
	}
	double Get_sig_pp_tot() const
	{
		return sig_pp_tot;
	}

	void Set_sig_pp_el(double a)
	{
		sig_pp_el = a;
	}
	double Get_sig_pp_el() const
	{
		return sig_pp_el;
	}

	void Set_sig_pn_el(double a)
	{
		sig_pn_el = a;
	}
	double Get_sig_pn_el() const
	{
		return sig_pn_el;
	}

	void Set_sig_pp_sd(double a)
	{
		sig_pp_sd = a;
	}
	double Get_sig_pp_sd() const
	{
		return sig_pp_sd;
	}

	void Set_sig_pn_sd(double a)
	{
		sig_pn_sd = a;
	}
	double Get_sig_pn_sd() const
	{
		return sig_pn_sd;
	}

	void Set_sig_pN_tot(double a)
	{
		sig_pN_tot = a;
	}
	double Get_sig_pN_tot() const
	{
		return sig_pN_tot;
	}

	void Set_sig_pN_inel(double a)
	{
		sig_pN_inel = a;
	}
	double Get_sig_pN_inel() const
	{
		return sig_pN_inel;
	}

	void Set_sig_pN_el(double a)
	{
		sig_pN_el = a;
	}
	double Get_sig_pN_el() const
	{
		return sig_pN_el;
	}

	void Set_lambda_tot(double a)
	{
		lambda_tot = a;
	}
	double Get_lambda_tot() const
	{
		return lambda_tot;
	}

	void Set_elastic_diff(double a)
	{
		elastic_diff = a;
	}
	double Get_elastic_diff() const
	{
		return elastic_diff;
	}

	void Set_com_sqd(double a)
	{
		com_sqd = a;
	}
	double Get_com_sqd() const
	{
		return com_sqd;
	}

	void Set_density(double a)
	{
		density = a;
	}
	double Get_density() const
	{
		return density;
	}

	void Set_atomic_mass(double a)
	{
		atomic_mass = a;
	}
	double Get_atomic_mass() const
	{
		return atomic_mass;
	}

	void Set_atomic_no(double a)
	{
		atomic_no = a;
	}
	double Get_atomic_no() const
	{
		return atomic_no;
	}

	void Set_scat_type(int a)
	{
		scat_type = a;
	}
	int Get_scat_type() const
	{
		return scat_type;
	}

	ParticleTracking::ppElasticScatter* GetElasticScatter() const
	{
		return ElasticScatter;
	}

	ParticleTracking::ppDiffractiveScatter* GetDiffractiveScatter() const
	{
		return DiffractiveScatter;
	}

private:
	ParticleTracking::ppElasticScatter* ElasticScatter;
	ParticleTracking::ppDiffractiveScatter* DiffractiveScatter;
	double E0;
	double sig_pN_tot_ref;
	double sig_pN_inel_ref;
	double sig_R_ref;
	double sig_R;
	double sig_pp_tot;
	double sig_pp_el;
	double sig_pn_el;
	double sig_pp_sd;
	double sig_pn_sd;
	double sig_pN_tot;
	double sig_pN_inel;
	double sig_pN_el;
	double lambda_tot;
	double elastic_diff;
	double com_sqd;
	double density;
	double atomic_mass;
	double atomic_no;
	int scat_type;
};
}

#endif
