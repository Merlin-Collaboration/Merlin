/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "AcceleratorModel.h"

#include "utils.h"

#include "RandomNG.h"

#include "CCFailureProcess.h"

#include "LatticeFunctions.h"
#include "PhaseAdvance.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace std;

namespace ParticleTracking
{

CCFailureProcess::CCFailureProcess(int priority, int mode, AcceleratorModel* model, LatticeFunctionTable* twiss) :
	ParticleBunchProcess("CRAB CAVITY FAILURE", priority), AccModelCC(model), TwissCC(twiss)
{
	ATLAS_on = 1;
	CMS_on = 1;
	Horizontal_CC = 0;
	upstream = 0;
	ATLAS = 0;
	IP1_up = 0;
	IP1_down = 0;
	IP5_up = 0;
	IP5_down = 0;
	IP1_up_count = 0;
	IP1_down_count = 0;
	IP5_up_count = 0;
	IP5_down_count = 0;
	Turn = 1;
	failure_on = 1;
	testn = 0;

	for(int i = 0; i < 4; i++)
	{
		Atlas_Upstream_V1[i] = 0.;
		CMS_Upstream_V1[i] = 0.;
		Atlas_Upstream_deltamu[i] = 0.;
		CMS_Upstream_deltamu[i] = 0.;
		Atlas_Upstream_elno[i] = 0.;
		CMS_Upstream_elno[i] = 0.;
		Atlas_Upstream_M12[i] = 0;
		CMS_Upstream_M12[i] = 0;
	}

	EnergyCC = 7E12;
	theta = 590e-6 / 2;           // Crossing angle
	omega = 400.79E6 * 2 * pi;  // Freq of CC
	phi_s = 0.0;                // CC phase - usually 0.0
	n = 4;
	fail_turns = 3;
	non_fail_turns = 10;
}

CCFailureProcess::CCFailureProcess(int priority, int mode, AcceleratorModel* model, LatticeFunctionTable* twiss, double
	freq, double crossing, double phase) :
	ParticleBunchProcess("CRAB CAVITY FAILURE", priority), AccModelCC(model), TwissCC(twiss), omega(freq), theta(
		crossing), phi_s(phase)
{
	ATLAS_on = 1;
	CMS_on = 1;
	Horizontal_CC = 0;
	upstream = 0;
	ATLAS = 0;
	IP1_up = 0;
	IP1_down = 0;
	IP5_up = 0;
	IP5_down = 0;
	IP1_up_count = 0;
	IP1_down_count = 0;
	IP5_up_count = 0;
	IP5_down_count = 0;
	Turn = 1;
	EnergyCC = 7E12;
	n = 4;
	failure_on = 1;
	testn = 0;

	for(int i = 0; i < 4; i++)
	{
		Atlas_Upstream_V1[i] = 0.;
		CMS_Upstream_V1[i] = 0.;
		Atlas_Upstream_deltamu[i] = 0.;
		CMS_Upstream_deltamu[i] = 0.;
		Atlas_Upstream_elno[i] = 0.;
		CMS_Upstream_elno[i] = 0.;
	}

	fail_turns = 3;
	non_fail_turns = 10;
}

CCFailureProcess::CCFailureProcess(int priority, int mode, AcceleratorModel* model, LatticeFunctionTable* twiss, double
	freq, double crossing, double phase, int non_fail_turn, int fail_turn) :
	ParticleBunchProcess("CRAB CAVITY FAILURE", priority), AccModelCC(model), TwissCC(twiss), omega(freq), theta(
		crossing), phi_s(phase), non_fail_turns(non_fail_turn), fail_turns(fail_turn)
{
	ATLAS_on = 1;
	CMS_on = 1;
	Horizontal_CC = 0;
	upstream = 0;
	ATLAS = 0;
	IP1_up = 0;
	IP1_down = 0;
	IP5_up = 0;
	IP5_down = 0;
	IP1_up_count = 0;
	IP1_down_count = 0;
	IP5_up_count = 0;
	IP5_down_count = 0;
	Turn = 1;
	EnergyCC = 7E12;
	n = 4;
	failure_on = 1;
	testn = 0;

	for(int i = 0; i < 4; i++)
	{
		Atlas_Upstream_V1[i] = 0.;
		CMS_Upstream_V1[i] = 0.;
		Atlas_Upstream_deltamu[i] = 0.;
		CMS_Upstream_deltamu[i] = 0.;
		Atlas_Upstream_elno[i] = 0.;
		CMS_Upstream_elno[i] = 0.;
	}
}

void CCFailureProcess::InitialiseProcess(Bunch& bunch)
{
	ParticleBunchProcess::InitialiseProcess(bunch);
	if(!currentBunch)
	{
		cout << "CCFailure warning: !currentBunch" << endl;
	}
}

void CCFailureProcess::SetCurrentComponent(AcceleratorComponent& component)
{
	CrabMarker* aCrabMarker = dynamic_cast<CrabMarker*>(&component);
	active = (currentBunch != nullptr) && (aCrabMarker);

	if(active)
	{
		currentComponent = &component;

		if(!ATLAS_on && !CMS_on)
		{
			cout << "CCFailureProcess: Warning: Both ATLAS and CMS switches are off. Failure Process does nothing"
				 << endl;
		}

		else if(ATLAS_on && CMS_on)
		{
			// This checks if we have made a full turn of the machine in terms of Crab Cavities
			if(IP1_up && IP1_down && IP5_up && IP5_down)
			{
				Turn++;
				cout << "CrabFailureProcess ++TURN, Turn = " << Turn << endl;
				//Reset turn
				IP1_up = 0;
				IP1_down = 0;
				IP5_up = 0;
				IP5_down = 0;
				//Reset CC counter
				IP5_up_count = 0;
				IP5_down_count = 0;
				IP1_up_count = 0;
				IP1_down_count = 0;
			}
		}
		else if(!ATLAS_on && CMS_on)
		{
			// This checks if we have made a full turn of the machine in terms of Crab Cavities
			if(IP5_up && IP5_down)
			{
				Turn++;
				cout << "CrabFailureProcess ++TURN, Turn = " << Turn << endl;
				//Reset turn
				IP5_up = 0;
				IP5_down = 0;
				//Reset CC counter
				IP5_up_count = 0;
				IP5_down_count = 0;
			}
		}

		else if(ATLAS_on && !CMS_on)
		{
			// This checks if we have made a full turn of the machine in terms of Crab Cavities
			if(IP1_up && IP1_down)
			{
				Turn++;
				cout << "CrabFailureProcess ++TURN, Turn = " << Turn << endl;
				//Reset turn
				IP1_up = 0;
				IP1_down = 0;
				//Reset CC counter
				IP1_up_count = 0;
				IP1_down_count = 0;
			}
		}

		//Here we check which set of CCs we are using, and count each so that we can know when we have traversed all 16 for beam1
		//Note that in order to use process this we must inject a Gaussian bunch immediately before a set of upstream CCs
		if((currentComponent->GetComponentLatticePosition() >= 10000) &&
			((currentComponent->GetComponentLatticePosition() <= 13200)))
		{
			upstream = 1;
			ATLAS = 0;
			Horizontal_CC = 1;
			IP5_up_count++;

			if(IP5_up_count == 4)
			{
				IP5_up = 1;
			}
			else
			{
				IP5_up = 0;
			}
		}
		else if((currentComponent->GetComponentLatticePosition() >= 13300) &&
			((currentComponent->GetComponentLatticePosition() <= 15000)))
		{
			upstream = 0;
			ATLAS = 0;
			Horizontal_CC = 1;
			IP5_down_count++;

			if(IP5_down_count == 4)
			{
				IP5_down = 1;
			}
			else
			{
				IP5_down = 0;
			}
		}
		else if(currentComponent->GetComponentLatticePosition() <= 200)
		{
			upstream = 0;
			ATLAS = 1;
			Horizontal_CC = 0;
			IP1_down_count++;

			if(IP1_down_count == 4)
			{
				IP1_down = 1;
			}
			else
			{
				IP1_down = 0;
			}
		}
		else if((currentComponent->GetComponentLatticePosition() >= 20000) &&
			((currentComponent->GetComponentLatticePosition() <= 30000)))
		{
			upstream = 1;
			ATLAS = 1;
			Horizontal_CC = 0;
			IP1_up_count++;

			if(IP1_up_count == 4)
			{
				IP1_up = 1;
			}
			else
			{
				IP1_up = 0;
			}
		}

		Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
		Beta_p = LorentzBeta(Gamma_p);
	}
	else
	{
		currentComponent = nullptr;
	}
}

double CCFailureProcess::GetMaxAllowedStepSize() const
{
	return currentComponent->GetLength();
}

void CCFailureProcess::DoProcess(double ds)
{

	if(!ATLAS_on && !CMS_on)
	{
	}

	else if(ATLAS_on && CMS_on)
	{
		if(Turn < (non_fail_turns + fail_turns))
		{
			ParticleBunch* newbunch = new ParticleBunch(currentBunch->GetReferenceMomentum(),
				currentBunch->GetTotalCharge() / currentBunch->size());
			newbunch->clear();
			newbunch->reserve(currentBunch->size());

			if(Beta_p == 0)
			{
				Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
				Beta_p = LorentzBeta(Gamma_p);
			}

			double DeltaMuX(0.), DeltaMuY(0.), MuX(0.), MuY(0.), M12(0.), M22(0.), V1(0.), V2(0.);
			int n1(0.), n2(0.);

			++testn;
			cout << "\n\t CCFAILUREPROCESS Called " << testn << " times " << endl;

			// These depend on whether we are pre or post IP
			n1 = AccModelCC->FindElementLatticePosition(currentComponent->GetName());
			if(ATLAS)
			{
				n2 = AccModelCC->FindElementLatticePosition("IP1.L1");
			}
			else
			{
				n2 = AccModelCC->FindElementLatticePosition("IP5") + 1; //+1 as phase is incorrect
			}

			//Calc Mu / DeltaMu
			if(upstream)
			{
				pair<double, double> DeltaMu = CalcDeltaMu(n1, n2);
				DeltaMuX = sqrt(pow(DeltaMu.first, 2)) * 2 * pi;
				DeltaMuY = sqrt(pow(DeltaMu.second, 2)) * 2 * pi;
			}
			else
			{
				pair<double, double> DeltaMu;
				if(ATLAS)       //Vertical
				{
					DeltaMu = CalcMu(n1);

					DeltaMuX = sqrt(pow(DeltaMu.first, 2)) * 2 * pi;
					DeltaMuY = sqrt(pow(DeltaMu.second, 2)) * 2 * pi + Atlas_Upstream_deltamu[IP1_down_count - 1];
				}
				else            //Horizontal
				{
					DeltaMu = CalcDeltaMu(n1, n2);

					DeltaMuX = sqrt(pow(DeltaMu.first, 2)) * 2 * pi + CMS_Upstream_deltamu[IP1_down_count - 1];
					DeltaMuY = sqrt(pow(DeltaMu.second, 2)) * 2 * pi;
				}
			}

			//Calc M12 & M22
			if(Horizontal_CC)  //CMS
			{
				if(upstream)
				{
					M12 = CalcM_12(n1, n2, Horizontal_CC);
				}
				else
				{
					M12 = CalcM_12(n2, n1, DeltaMuX / 2, Horizontal_CC);
					M22 = CalcM_22(CMS_Upstream_elno[IP5_down_count - 1], n1, DeltaMuX, Horizontal_CC);
				}
			}
			else  //ATLAS
			{
				if(upstream)
				{
					M12 = CalcM_12(n1, n2, Horizontal_CC);
				}
				else
				{
					M12 = CalcM_12(n2, n1, DeltaMuY / 2, Horizontal_CC);
					M22 = CalcM_22(Atlas_Upstream_elno[IP1_down_count - 1], n1, DeltaMuY, Horizontal_CC);
				}
			}

			//Store V1s
			if(ATLAS && upstream)
			{
				V1 = CalcV1(M12);

				Atlas_Upstream_elno[IP1_up_count - 1] = n1;
				Atlas_Upstream_V1[IP1_up_count - 1] = V1;
				Atlas_Upstream_deltamu[IP1_up_count - 1] = DeltaMuY;
				Atlas_Upstream_M12[IP1_up_count - 1] = M12;
			}
			else if(upstream)
			{
				V1 = CalcV1(M12);

				CMS_Upstream_elno[IP5_up_count - 1] = n1;
				CMS_Upstream_V1[IP5_up_count - 1] = V1;
				CMS_Upstream_deltamu[IP5_up_count - 1] = DeltaMuX;
				CMS_Upstream_M12[IP5_up_count - 1] = M12;
			}
			//Load V1s
			if(ATLAS && !upstream)
			{
				V2 = CalcV2(Atlas_Upstream_V1[IP1_down_count - 1], M22);
			}
			else if(!upstream)
			{
				V2 = CalcV2(CMS_Upstream_V1[IP5_down_count - 1], M22);
			}

//FAILURE
			double fail_interval = 1 / (double) fail_turns;
			if(upstream && ATLAS && (Turn >= non_fail_turns) && (failure_on))
			{
				V1 = V1 * (1 - (((Turn + 1) - non_fail_turns) * fail_interval));
			}

			cout << "\nElement = " << currentComponent->GetName() << ", position = "
				 << currentComponent->GetComponentLatticePosition() << endl;
			if(upstream)
			{
				if(ATLAS)
				{
					cout << "Upstream ATLAS V1 = " << V1 << " \nDeltaMuY = " << DeltaMuY * (360 / (2 * pi))
						 << " M12 = " << M12 << endl;
				}
				else
				{
					cout << "Upstream CMS V1 = " << V1 << " \nDeltaMuX = " << DeltaMuX * (360 / (2 * pi))
						 << " M12 = " << M12 << endl;
				}
			}
			else
			{
				if(ATLAS && (Atlas_Upstream_V1[3] != 0.))
				{
					cout << "Downstream ATLAS V2 = " << V2 << "\nDeltaMuY = " << DeltaMuY * (360 / (2 * pi))
						 << " M12 = " << M12 << " M22 = " << M22 << endl;
				}
				else if(!ATLAS && (CMS_Upstream_V1[3] != 0.))
				{
					cout << "Downstream CMS V2 = " << V2 << "\nDeltaMuX = " << DeltaMuX * (360 / (2 * pi))
						 << " M12 = " << M12 << " M22 = " << M22 << endl;
				}

			}

			for(PSvectorArray::iterator p = currentBunch->begin(); p != currentBunch->end(); p++)
			{
				if(ATLAS)
				{
					if(upstream)
					{
						ApplyPreCCKick(*p, V1, M12, Horizontal_CC);
					}
					else
					{
						ApplyPostCCKick(*p, V2, M12, Horizontal_CC);
					}
				}
				if(std::isnan(p->x()) || std::isnan(p->xp()) || std::isnan(p->y()) || std::isnan(p->yp()) || std::isnan(
						p->ct()) || std::isnan(p->dp()))
				{
					cout << "CC Particle Lost" << endl;
				}
				else
				{
					newbunch->AddParticle(*p);
				}
			}

			//Apply changes to particle bunch
			currentBunch->clear();
			currentBunch->swap(*newbunch);
		}
	}

	else if(ATLAS_on && !CMS_on && ATLAS)
	{
		if(Turn < (non_fail_turns + fail_turns))
		{
			ParticleBunch* newbunch = new ParticleBunch(currentBunch->GetReferenceMomentum(),
				currentBunch->GetTotalCharge() / currentBunch->size());
			newbunch->clear();
			newbunch->reserve(currentBunch->size());

			if(Beta_p == 0)
			{
				Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
				Beta_p = LorentzBeta(Gamma_p);
			}

			double DeltaMuY(0.), MuX(0.), MuY(0.), M12(0.), M22(0.), V1(0.), V2(0.);
			int n1(0.), n2(0.);

			++testn;
			cout << "\n\t CCFAILUREPROCESS Called " << testn << " times " << endl;

			// These depend on whether we are pre or post IP
			n1 = AccModelCC->FindElementLatticePosition(currentComponent->GetName());
			n2 = AccModelCC->FindElementLatticePosition("IP1.L1");

			//Calc Mu / DeltaMu
			if(upstream)
			{
				pair<double, double> DeltaMu = CalcDeltaMu(n1, n2);
				DeltaMuY = sqrt(pow(DeltaMu.second, 2)) * 2 * pi;
			}
			else
			{
				pair<double, double> DeltaMu;
				DeltaMu = CalcMu(n1);
				DeltaMuY = sqrt(pow(DeltaMu.second, 2)) * 2 * pi + Atlas_Upstream_deltamu[IP1_down_count - 1];
			}

			//Calc M12 & M22
			if(upstream)
			{
				M12 = CalcM_12(n1, n2, Horizontal_CC);
			}
			else
			{
				M12 = CalcM_12(n2, n1, DeltaMuY / 2, Horizontal_CC);
				M22 = CalcM_22(Atlas_Upstream_elno[IP1_down_count - 1], n1, DeltaMuY, Horizontal_CC);
			}

			//Store V1s
			if(upstream)
			{
				V1 = CalcV1(M12);

				Atlas_Upstream_elno[IP1_up_count - 1] = n1;
				Atlas_Upstream_V1[IP1_up_count - 1] = V1;
				Atlas_Upstream_deltamu[IP1_up_count - 1] = DeltaMuY;
				Atlas_Upstream_M12[IP1_up_count - 1] = M12;
			}
			//Load V1s
			else if(!upstream)
			{
				V2 = CalcV2(Atlas_Upstream_V1[IP1_down_count - 1], M22);
			}

//FAILURE
			double fail_interval = 1 / (double) fail_turns;
			if(!upstream && (Turn >= non_fail_turns) && (failure_on))
			{
				//140 degrees is 2.443 radians 140 * (2*pi)/360
				double phi = 90 * (2 * pi) / 360;
				phi_s = (((Turn) - (non_fail_turns - 1)) * fail_interval) * phi;
			}

			cout << "\nElement = " << currentComponent->GetName() << ", position = "
				 << currentComponent->GetComponentLatticePosition() << endl;
			if(upstream)
			{
				if(ATLAS)
				{
					cout << "Upstream ATLAS V1 = " << V1 << " \nDeltaMuY = " << DeltaMuY * (360 / (2 * pi))
						 << " M12 = " << M12 << " Phi = " << phi_s << endl;
				}
			}
			else
			{
				if(ATLAS && (Atlas_Upstream_V1[3] != 0.))
				{
					cout << "Downstream ATLAS V2 = " << V2 << "\nDeltaMuY = " << DeltaMuY * (360 / (2 * pi))
						 << " M12 = " << M12 << " M22 = " << M22 << " Phi = " << phi_s << endl;
				}

			}

			for(PSvectorArray::iterator p = currentBunch->begin(); p != currentBunch->end(); p++)
			{
				if(upstream)
				{
					ApplyPreCCKick(*p, V1, M12, Horizontal_CC);
				}
				else
				{
					ApplyPostCCKick(*p, V2, M12, Horizontal_CC);
				}

				if(std::isnan(p->x()) || std::isnan(p->xp()) || std::isnan(p->y()) || std::isnan(p->yp()) || std::isnan(
						p->ct()) || std::isnan(p->dp()))
				{
					cout << "CC Particle Lost" << endl;
				}
				else
				{
					newbunch->AddParticle(*p);
				}
			}

			//Apply changes to particle bunch
			currentBunch->clear();
			currentBunch->swap(*newbunch);
			phi_s = 0.0;
		}
	}

	else if(!ATLAS_on && CMS_on && !ATLAS)
	{
		if(Turn < (non_fail_turns + fail_turns))
		{
			ParticleBunch* newbunch = new ParticleBunch(currentBunch->GetReferenceMomentum(),
				currentBunch->GetTotalCharge() / currentBunch->size());
			newbunch->clear();
			newbunch->reserve(currentBunch->size());

			if(Beta_p == 0)
			{
				Gamma_p = LorentzGamma(currentBunch->GetReferenceMomentum(), ProtonMass);
				Beta_p = LorentzBeta(Gamma_p);
			}

			double DeltaMuX(0.), MuX(0.), MuY(0.), M12(0.), M22(0.), V1(0.), V2(0.);

			int n1(0.), n2(0.);

			++testn;
			cout << "\n\t CCFAILUREPROCESS Called " << testn << " times " << endl;

			// These depend on whether we are pre or post IP
			n1 = AccModelCC->FindElementLatticePosition(currentComponent->GetName());
			n2 = AccModelCC->FindElementLatticePosition("IP5") + 1; //+1 as phase is incorrect

			//Calc Mu / DeltaMu

			pair<double, double> DeltaMu = CalcDeltaMu(n1, n2);
			if(upstream)
			{
				DeltaMuX = sqrt(pow(DeltaMu.first, 2)) * 2 * pi;
			}
			else
			{
				DeltaMuX = sqrt(pow(DeltaMu.first, 2)) * 2 * pi + CMS_Upstream_deltamu[IP5_down_count - 1];
			}

			//Calc M12 & M22
			if(upstream)
			{
				M12 = CalcM_12(n1, n2, Horizontal_CC);
			}
			else
			{
				M12 = CalcM_12(n2, n1, DeltaMuX / 2, Horizontal_CC);
				M22 = CalcM_22(CMS_Upstream_elno[IP5_down_count - 1], n1, DeltaMuX, Horizontal_CC);
			}

			//Store V1s
			if(upstream)
			{
				V1 = CalcV1(M12);

				CMS_Upstream_elno[IP5_up_count - 1] = n1;
				CMS_Upstream_V1[IP5_up_count - 1] = V1;
				CMS_Upstream_deltamu[IP5_up_count - 1] = DeltaMuX;
				CMS_Upstream_M12[IP5_up_count - 1] = M12;
			}
			//Load V1s
			else
			{
				V2 = CalcV2(CMS_Upstream_V1[IP5_down_count - 1], M22);
			}

//FAILURE
			double fail_interval = 1 / (double) fail_turns;
			if(upstream && (Turn >= non_fail_turns) && (failure_on))
			{
				V1 = V1 * (1 - (((Turn + 1) - non_fail_turns) * fail_interval));
			}

			cout << "\nElement = " << currentComponent->GetName() << ", position = "
				 << currentComponent->GetComponentLatticePosition() << endl;
			if(upstream)
			{
				cout << "Upstream CMS V1 = " << V1 << " \nDeltaMuX = " << DeltaMuX * (360 / (2 * pi)) << " M12 = "
					 << M12 << endl;
			}
			else
			{
				cout << "Downstream CMS V2 = " << V2 << "\nDeltaMuX = " << DeltaMuX * (360 / (2 * pi)) << " M12 = "
					 << M12 << " M22 = " << M22 << endl;
			}

			for(PSvectorArray::iterator p = currentBunch->begin(); p != currentBunch->end(); p++)
			{
				if(upstream)
				{
					ApplyPreCCKick(*p, V1, M12, Horizontal_CC);
				}
				else
				{
					ApplyPostCCKick(*p, V2, M12, Horizontal_CC);
				}

				if(std::isnan(p->x()) || std::isnan(p->xp()) || std::isnan(p->y()) || std::isnan(p->yp()) || std::isnan(
						p->ct()) || std::isnan(p->dp()))
				{
					cout << "CC Particle Lost" << endl;
				}
				else
				{
					newbunch->AddParticle(*p);
				}
			}

			//Apply changes to particle bunch
			currentBunch->clear();
			currentBunch->swap(*newbunch);
		}
	}
}

double CCFailureProcess::CalcM_12(int start, int end, double deltamu, bool horizontal)
{
	double beta1(0.), beta2(0.);

	if(horizontal)
	{
		beta1 = TwissCC->Value(1, 1, 1, start);
		beta2 = TwissCC->Value(1, 1, 1, end);
	}
	else
	{
		beta1 = TwissCC->Value(3, 3, 2, start);
		beta2 = TwissCC->Value(3, 3, 2, end);
	}

	return sqrt(beta1 * beta2) * sin(deltamu);
}

double CCFailureProcess::CalcM_22(int start, int end, double deltamu, bool horizontal)
{
	//Simplified calculation omitting Alpha_2*sin(deltamu) term
	double beta1(0.), beta2(0.);

	if(horizontal)
	{
		beta1 = TwissCC->Value(1, 1, 1, start);
		beta2 = TwissCC->Value(1, 1, 1, end);
	}
	else
	{
		beta1 = TwissCC->Value(3, 3, 2, start);
		beta2 = TwissCC->Value(3, 3, 2, end);
	}

	return sqrt(beta1 / beta2) * cos(deltamu);
}

double CCFailureProcess::CalcM_12(int start, int end, bool horizontal)
{
	double beta1(0.), beta2(0.), deltamu(0.);
	pair<double, double> DeltaMu = CalcDeltaMu(start, end);
	double DeltaMuX = sqrt(pow(DeltaMu.first, 2)) * 2 * pi;
	double DeltaMuY = sqrt(pow(DeltaMu.second, 2)) * 2 * pi;

	if(horizontal)
	{
		deltamu = DeltaMuX;
		beta1 = TwissCC->Value(1, 1, 1, start);
		beta2 = TwissCC->Value(1, 1, 1, end);
	}
	else
	{
		deltamu = DeltaMuY;
		beta1 = TwissCC->Value(3, 3, 2, start);
		beta2 = TwissCC->Value(3, 3, 2, end);
	}

	return sqrt(beta1 * beta2) * sin(deltamu);
}

double CCFailureProcess::CalcM_22(int start, int end, bool horizontal)
{
	//Simplified calculation omitting Alpha_2*sin(deltamu) term
	double beta1(0.), beta2(0.), deltamu(0.);
	pair<double, double> DeltaMu = CalcDeltaMu(start, end);
	double DeltaMuX = sqrt(pow(DeltaMu.first, 2)) * 2 * pi;
	double DeltaMuY = sqrt(pow(DeltaMu.second, 2)) * 2 * pi;

	if(horizontal)
	{
		deltamu = DeltaMuX;
		beta1 = TwissCC->Value(1, 1, 1, start);
		beta2 = TwissCC->Value(1, 1, 1, end);
	}
	else
	{
		deltamu = DeltaMuY;
		beta1 = TwissCC->Value(3, 3, 2, start);
		beta2 = TwissCC->Value(3, 3, 2, end);
	}

	return sqrt(beta1 / beta2) * cos(deltamu);
}

pair<double, double> CCFailureProcess::CalcMu(int element)
{
	PhaseAdvance* PA = new PhaseAdvance(AccModelCC, TwissCC, EnergyCC);
	pair<double, double> Mu = PA->CalcIntegerPart(element);
	return Mu;
}

pair<double, double> CCFailureProcess::CalcDeltaMu(int element1, int element2)
{
	pair<double, double> Mu1 = CalcMu(element1);
	pair<double, double> Mu2 = CalcMu(element2);

	double delta_mu_x = Mu2.first - Mu1.first;
	double delta_mu_y = Mu2.second - Mu1.second;

	pair<double, double> DeltaMu = std::make_pair(delta_mu_x, delta_mu_y);

	return DeltaMu;
}

double CCFailureProcess::CalcV1(double M12)
{
	return (SpeedOfLight * EnergyCC * tan(theta)) / (omega * M12 * n);
}

double CCFailureProcess::CalcV1(double deltamu, int n1, int n2, bool horizontal)
{
	if(horizontal)
	{
		return (SpeedOfLight * EnergyCC * tan(theta) * 1E-6) / (omega * sqrt(TwissCC->Value(1, 1, 1, n1)
			   * TwissCC->Value(1, 1, 1, n2)) * sin(deltamu));
	}
	else
	{
		return (SpeedOfLight * EnergyCC * tan(theta) * 1E-6) / (omega * sqrt(TwissCC->Value(3, 3, 2, n1)
			   * TwissCC->Value(3, 3, 2, n2)) * sin(deltamu));
	}
}

double CCFailureProcess::CalcV2(double V1, double M22)
{
	return -1 * M22 * V1;
}

void CCFailureProcess::ApplyPreCCKick(PSvector &p, double V, double M12, bool horizontal)
{
	double omega_ov_c = omega / SpeedOfLight;
	double tantheta = tan(theta);
	int n = 4;  //No of CCs

	//Working
	double kick = (V * sin(phi_s + (p.ct() * omega_ov_c))) / (EnergyCC);

	if(horizontal)
	{
		p.xp() -= kick;
	}
	else
	{
		p.yp() -= kick;
	}
}

void CCFailureProcess::ApplyPostCCKick(PSvector &p, double V, double M12, bool horizontal)
{
	double omega_ov_c = omega / SpeedOfLight;
	double tantheta = tan(theta);
	int n = 4;  //No of CCs

	//Working
	double kick = (V * sin(phi_s + (p.ct() * omega_ov_c))) / (EnergyCC);

	if(horizontal)
	{
		p.xp() -= kick;
	}
	else
	{
		p.yp() -= kick;
	}
}

} //End namespace ParticleTracking
