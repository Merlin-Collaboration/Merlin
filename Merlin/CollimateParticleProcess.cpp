/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iterator>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <sstream>

#include "merlin_config.h"

#include "Aperture.h"
#include "InterpolatedApertures.h"
#include "CollimatorAperture.h"
#include "Collimator.h"

#include "ParticleComponentTracker.h"

#include "CollimateParticleProcess.h"

#include "utils.h"
#include "PhysicalUnits.h"

namespace
{

using namespace ParticleTracking;

void OutputIndexParticles(const PSvectorArray lost_p, const list<size_t>& lost_i, ostream& os)
{
	PSvectorArray::const_iterator p = lost_p.begin();
	list<size_t>::const_iterator ip = lost_i.begin();

	while(p != lost_p.end())
	{
		os << std::setw(12) << right << *ip;
		os << *p;
		++p;
		++ip;
	}
}

} // end anonymous namespace

namespace ParticleTracking
{

CollimateParticleProcess::CollimateParticleProcess(int priority, int mode, std::ostream* osp) :
	ParticleBunchProcess("PARTICLE COLLIMATION", priority), cmode(mode), os(osp), createLossFiles(false), file_prefix(
		""), lossThreshold(1), nstart(0), pindex(nullptr), CollimationOutputSet(false), ColParProTurn(0),
	FirstElementSet(0), scatter(false), bin_size(0.1 * PhysicalUnits::meter), Imperfections(false)
{
}

CollimateParticleProcess::~CollimateParticleProcess()
{
	if(pindex != nullptr)
	{
		delete pindex;
	}
}

void CollimateParticleProcess::InitialiseProcess(Bunch& bunch)
{
	ParticleBunchProcess::InitialiseProcess(bunch);
	idtbl.clear();
	if(currentBunch)
	{
		nstart = currentBunch->size();
		nlost = 0;
		if(pindex != nullptr)
		{
			pindex->clear();
			for(size_t n = 0; n < nstart; n++)
			{
				pindex->push_back(n);
			}
		}
	}
}

void CollimateParticleProcess::SetCurrentComponent(AcceleratorComponent& component)
{
	if(!FirstElementSet)
	{
		FirstElementName = component.GetName();
		FirstElementS = component.GetComponentLatticePosition();
		FirstElementSet = 1;
		ColParProTurn = 1;
	}
	else if(component.GetName() == FirstElementName && component.GetComponentLatticePosition() == FirstElementS)
	{
		++ColParProTurn;
	}

	active = (currentBunch != nullptr) && (component.GetAperture() != nullptr);
	if(active)
	{
		currentComponent = &component;
		s = 0;
		Collimator* aCollimator = dynamic_cast<Collimator*>(&component);

		const CollimatorAperture* tap = dynamic_cast<const CollimatorAperture*> (currentComponent->GetAperture());
		is_collimator = scatter && tap;

		if(!is_collimator)
		{
			// not a collimator, so set up for normal hard-edge collimation
			at_entr = (COLL_AT_ENTRANCE & cmode) != 0;
			at_cent = (COLL_AT_CENTER & cmode) != 0;
			at_exit = (COLL_AT_EXIT & cmode) != 0;
			SetNextS();
		}
		else
		{
			at_cent = at_entr = false; // currently scatter only at exit
			at_exit = true;
			SetNextS();
			currentBunch->SetScatterConfigured(false);
			len = aCollimator->GetLength();
			//	CollimatorAperture* CollimatorJaw = dynamic_cast<CollimatorAperture*>(aCollimator->GetAperture());
		}

		//For precision tracking of lost particles in non-collimators
		//make a copy of the input array if we are a magnet.
		//This should also occur before any tracking, hence also any poleface rotations on a dipole, etc
		if(!is_collimator)
		{
			InputArray = currentBunch->GetParticles();
		}
	}
	else
	{
		s_total += component.GetLength();
		currentComponent = nullptr;
	}
}

void CollimateParticleProcess::DoProcess(double ds)
{
	s += ds;

	if(fequal(s, next_s))
	{
		//This lets the scattering routine know how far down the collimator we are for aperture checking inside the scattering step.
		currentBunch->SetIntS(s - ds);
		DoCollimation();
		SetNextS();
	}

	// If we are finished, GetNextS() will have set the process inactive.
	// In that case we can update s_total with the component length.
	if(!active)
	{
		s_total += currentComponent->GetLength();
	}
}

double CollimateParticleProcess::GetMaxAllowedStepSize() const
{
	if(!is_collimator)
	{
		return next_s - s;
	}
	else
	{
		return bin_size;
	}
}

void CollimateParticleProcess::IndexParticles(bool index)
{
	if(index && pindex == nullptr)
	{
		pindex = new list<size_t>;
	}
	else if(!index && pindex != nullptr)
	{
		delete pindex;
		pindex = nullptr;
	}
}

void CollimateParticleProcess::IndexParticles(list<size_t>& anIndex)
{
	if(!pindex)
	{
		delete pindex;
	}

	pindex = &anIndex;
}

void CollimateParticleProcess::SetLossThreshold(double losspc)
{
	lossThreshold = losspc / 100.0;
}

void CollimateParticleProcess::DoCollimation()
{
	//The aperture of this element
	const Aperture *ap = currentComponent->GetAperture();

	// If there are no losses there is no need to go through the expensive
	// process of copying all the particles to a new bunch. So check first
	bool any_loss = false;
	size_t first_loss = 0;
	if(is_collimator)
	{
		for(PSvectorArray::iterator p = currentBunch->begin(); p != currentBunch->end(); p++)
		{
			if(!ap->CheckWithinApertureBoundaries((*p).x() - bin_size * (*p).xp(), (*p).y() - bin_size * (*p).yp(), s))
			{
				any_loss = true;
				first_loss = p - currentBunch->begin();
				break;
			}
		}
	}
	else
	{
		for(PSvectorArray::iterator p = currentBunch->begin(); p != currentBunch->end(); p++)
		{
			if(!ap->CheckWithinApertureBoundaries((*p).x(), (*p).y(), s))
			{
				any_loss = true;
				first_loss = p - currentBunch->begin();
				break;
			}
		}
	}

	if(!any_loss)
	{
		return;
	}

	//The array of lost particles
	PSvectorArray lost;
	list<size_t> lost_i;

	list<size_t>::iterator ip;
	if(pindex != nullptr)
	{
		ip = pindex->begin();
	}

	//For copying surviving particles to, which is faster than deleting the individual lost particles
	ParticleBunch* NewBunch = new ParticleBunch(currentBunch->GetReferenceMomentum(), currentBunch->GetTotalCharge()
		/ currentBunch->size());
	NewBunch->reserve(currentBunch->size());

	size_t particle_number = 0;

	if(is_collimator)
	{
		for(PSvectorArray::iterator p = currentBunch->begin(); p != currentBunch->end();)
		{
			(*p).x() -= bin_size * (*p).xp();
			(*p).y() -= bin_size * (*p).yp();
			p++;
		}
	}

	for(PSvectorArray::iterator p = currentBunch->begin(); p != currentBunch->end();)
	{
		// If we are collimating at the end of the element, track back a drift
		// Do not do this at the start of the element.

//		if(is_collimator)
//		{
//			(*p).x() -= bin_size * (*p).xp();
//			(*p).y() -= bin_size * (*p).yp();
//		}
		if(particle_number >= first_loss && !ap->CheckWithinApertureBoundaries((*p).x(), (*p).y(), s))
		{
			// If the 'aperture' is a collimator, then the particle is lost
			// if the DoScatter(*p) returns true (energy cut)
			// If not a collimator, then do not scatter and directly remove the particle.
			if(!is_collimator || DoScatter(*p))
			{
				if(is_collimator)
				{
					(*p).ct() += (s - bin_size);
				}

				lost.push_back(*p);

				// This is slow for a STL Vector - instead we place the surviving particles into a new bunch and then swap - this is faster
				p++;
				if(pindex != nullptr)
				{
					lost_i.push_back(*ip);
					ip = pindex->erase(ip);
				}

				LostParticlePositions.push_back(particle_number);
			}
			else
			{
				(*p).location() = currentComponent->GetComponentLatticePosition();
				//Particle survives collimator
				NewBunch->AddParticle(*p);
				// need to increment iterators
				p++;
				if(pindex != nullptr)
				{
					ip++;
				}
			}
		}
		else
		{
			if(is_collimator)
			{
				(*p).x() += bin_size * (*p).xp();
				(*p).y() += bin_size * (*p).yp();
			}

			//Not interacting with the collimator: "Inside" the aperture; particle lives
			NewBunch->AddParticle(*p);
			p++;
			if(pindex != nullptr)
			{
				ip++;
			}
		}
		particle_number++;
	}

	currentBunch->clear();
	currentBunch->swap(*NewBunch);
	delete NewBunch;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Only copy the output if we are not a collimator and there are lost particles
	if(LostParticlePositions.size() != 0 && !is_collimator)
	{
		//make a new particle bunch to track the lost particles
		ParticleBunch* LostBunch = new ParticleBunch(currentBunch->GetReferenceMomentum(),
			currentBunch->GetTotalCharge() / currentBunch->size());
		double length = currentComponent->GetLength();
		//If we are dealing with a non-zero length element, we must do tracking
		if(length != 0)
		{
			//Clear out the old lost particles, these will all be at the end of the element which we do not want for a magnet.
			lost.clear();

			//Grab the lost particles from the copied input particle array and add them to the new particle bunch
			for(vector<unsigned int>::iterator p = LostParticlePositions.begin(); p != LostParticlePositions.end(); p++)
			{
				LostBunch->AddParticle(InputArray[*p]);
			}
			//Create a new tracker
			ParticleComponentTracker* LostParticleTracker = new ParticleComponentTracker();

			//Tell the new tracker to use the particle bunch of lost particles that has been made
			LostParticleTracker->SetBunch(*LostBunch);

			//Prepare the tracker to use the current accelerator component
			currentComponent->PrepareTracker(*LostParticleTracker);

			//While we are still inside the component
			while((LostParticleTracker->GetRemainingLength()) >= 0 && LostBunch->size() != 0)
			{
				double StepSize = bin_size;
				//If the remaining length of component is less than the step size, set the step size to this value
				if((LostParticleTracker->GetRemainingLength() - StepSize) < 0)
				{
					StepSize = LostParticleTracker->GetRemainingLength();
				}

				//Track the appropriate step length

				double IntegratedLength = LostParticleTracker->GetIntegratedLength();
				//Now loop over each particle in turn
				for(PSvectorArray::iterator p = LostBunch->begin(); p != LostBunch->end();)
				{
					//Check if the particle is outside the aperture
					//s, is where the integrator will start
					//LostParticleTracker->GetIntegratedLength() will give the position integrated past this point
					//(*p).ct() will give the offset for this specific particle

					if(!ap->CheckWithinApertureBoundaries((*p).x(), (*p).y(), IntegratedLength))
					{
						//if not, delete the particle, and add the coordinates to the lost bunch list (PSvectorArray lost)

						//Also set p.ct() as the length along the element!
						(*p).ct() += IntegratedLength;
						if((*p).ct() < 0)
						{
							(*p).ct() = 0;
						}
						if((*p).ct() > length)
						{

							(*p).ct() = length;
						}

						lost.push_back(*p);

						//CollimationOutput loss
						if(CollimationOutputSet && !is_collimator)
						{
							for(CollimationOutputIterator = CollimationOutputVector.begin();
								CollimationOutputIterator != CollimationOutputVector.end();
								++CollimationOutputIterator)
							{
								(*CollimationOutputIterator)->Dispose(*currentComponent, IntegratedLength, (*p),
									ColParProTurn);
							}
						}
						p = LostBunch->erase(p);
					}
					//else, the particle is inside and can be kept for this step
					else
					{
						p++;
					}
				}

				//Now move forward...
				if((LostParticleTracker->GetRemainingLength()) > 0)
				{
					LostParticleTracker->TrackStep(StepSize);
				}
				if(LostParticleTracker->GetRemainingLength() == 0)
				{
					break;
				}
			}

			//If there is anything left - possible bug.
			if(LostBunch->size() != 0)
			{
				std::cout << "POSSIBLE BUG: Leftovers: " << LostBunch->size() << "\t"
						  << currentComponent->GetQualifiedName() << "\t"
						  << LostParticleTracker->GetIntegratedLength() << "\t" << length << std::endl;
				for(PSvectorArray::iterator p = LostBunch->begin(); p != LostBunch->end(); p++)
				{
					(*p).ct() += LostParticleTracker->GetIntegratedLength();
					if((*p).ct() > length)
					{
						(*p).ct() = length;
					}
					lost.push_back(*p);
				}
			}

			//clean up the tracker
			delete LostParticleTracker;
			//and the "lost particle bunch"
			delete LostBunch;
		}
		//if the element has zero length nothing needs to be done since all the losses will have occurred at the same point anyway.
		//So PSvectorArray loss will contain the correct information
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	nlost += lost.size();
	// Old loss output - depreciated due to CollimationOutput
	//DoOutput(lost,lost_i);

	//make sure to clear up
	InputArray.clear();             //The input array
	LostParticlePositions.clear();  //A list of particles we want to use in the input array

	if(double(nlost) / double(nstart) >= lossThreshold)
	{
		std::cout << "nlost: " << nlost << "\tnstart: " << nstart << std::endl;
		throw ExcessiveParticleLoss(currentComponent->GetQualifiedName(), lossThreshold, nlost, nstart);
	}
}

void CollimateParticleProcess::SetNextS()
{

	if(at_entr)
	{
		next_s = 0;
		at_entr = false;
	}
	else if(at_cent)
	{
		next_s = currentComponent->GetLength() / 2;
		at_cent = false;
	}
	else if(at_exit)
	{
		next_s = currentComponent->GetLength();
		at_exit = false;
	}
	else
	{
		active = false;
	}

	if(is_collimator)
	{
		active = true;
		if((s + bin_size) > currentComponent->GetLength())
		{
			next_s = currentComponent->GetLength();
		}
		else
		{
			next_s = s + bin_size;
		}
	}
}

void CollimateParticleProcess::DoOutput(const PSvectorArray& lostb, const list<size_t>& lost_i)
{

	// Create a file and dump the lost particles
	// (if there are any)
	if(!lostb.empty())
	{
//	PSvectorArray lostp = bin_lost_output(lostb);
		if(os != nullptr)
		{
			double length = currentComponent->GetLength();
			double** lostp;
			int n;
			//deal with zero length elements separately
			if(length != 0)
			{
				//We bin in 0.1m sections of the element, this is how many bins we have
				//n is our number of bins
				n = (length / bin_size) + 1;

				bool overflow = false;

				//The element length may not be an exact multiple of the bin size, so we must take this into account.
				if((length / (double) n) != 0.0)
				{
					//Add an additional bin
					n++;
					overflow = true;
				}
				//Create the array - n rows, 3 cols
				lostp = new double*[n];
				for(int i = 0; i < n; ++i)
				{
					lostp[i] = new double[3];
				}

				//Initialize the array elements
				for(int j = 0; j < n; j++)
				{
					lostp[j][0] = bin_size * j;   //Start position
					lostp[j][1] = bin_size;     //Length
					lostp[j][2] = 0.0;      //Entries
				}

				//deal with the last bin if needed - will have a different length.
				if(overflow)
				{
					lostp[n - 1][0] = bin_size * (n - 1);
					lostp[n - 1][1] = length - (bin_size * (n - 1));
					lostp[n - 1][2] = 0.0;
				}

				for(size_t l = 0; l < lostb.size(); l++)
				{
					int x = (lostb[l].ct()) / bin_size;
					//Add one loss count to this bin
					lostp[x][2]++;
				}
			}
			//now deal with zero length elements;
			else
			{
				n = 1;
				lostp = new double*[1];
				lostp[0] = new double[3];
				lostp[0][0] = 0.0;      //Start position
				lostp[0][1] = bin_size;     //Length
				lostp[0][2] = lostb.size(); //Entries
			}

			//Now to do the output - first loop over each bin
			for(int j = 0; j < n; j++)
			{
				//We then check if there are any lost particles in this bin
				//If there are lost particles, then we must output the count (lost count > 0)
				if(lostp[j][2] > 0.0)
				{
					(*os).precision(16);
					(*os) << std::setw(35) << left << (*currentComponent).GetQualifiedName().c_str();           //Component name - can be quite long for certain LHC magnets
					//(*os) << std::setw(24)<<left<<lostp[j][0] + currentBunch->GetReferenceTime()-length;	//Bin start position
					(*os) << std::setw(24) << left << lostp[j][0] + currentComponent->GetComponentLatticePosition();    //Bin start position
					(*os) << std::setw(24) << left << lostp[j][1];                          //Bin length
					(*os) << "\t" << lostp[j][2];                                   //Loss count
					(*os) << "\t" << j * lostp[0][1];                               //Loss position in element
					(*os) << std::endl;
				}
			}
			delete[] lostp;
		}
		if(createLossFiles)
		{
			std::string id = (*currentComponent).GetQualifiedName();
			pair<IDTBL::iterator, bool> result = idtbl.insert(IDTBL::value_type(id, 0));
			int n = ++(*(result.first)).second;
			ostringstream fname;

			if(!file_prefix.empty())
			{
				fname << file_prefix;
			}

			fname << "_" << currentComponent->GetComponentLatticePosition() << "_";
			fname << id << '.' << n << ".loss";

			std::ofstream file(fname.str().c_str());

			if(!file)
			{
				std::cerr << "CollimateParticleProcess::DoOutput(): Failed to open " << fname.str() << std::endl;
				exit(EXIT_FAILURE);
			}

			if(pindex == nullptr)
			{
				copy(lostb.begin(), lostb.end(), ostream_iterator<PSvector>(file));
			}
			else
			{
				OutputIndexParticles(lostb, lost_i, file);
			}
		}
	}
}

bool CollimateParticleProcess::DoScatter(Particle& p)
{
	const CollimatorAperture *tap = (CollimatorAperture *) currentComponent->GetAperture();

	//int scatter_type = currentBunch->Scatter(p,len,tap);
	int scatter_type = currentBunch->Scatter(p, bin_size, tap);

	if(scatter_type == 1)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

ExcessiveParticleLoss::ExcessiveParticleLoss(const string& c_id, double threshold, size_t nlost, size_t nstart) :
	MerlinException()
{
	ostringstream buffer;
	buffer << "CollimateParticleProcess Exception\n";
	buffer << "particle loss threshold of " << 100 * threshold << "% exceeded ";
	buffer << '(' << nlost << '/' << nstart << ") at " << c_id;
	SetMsg(buffer.str());
}

void CollimateParticleProcess::SetOutputBinSize(double binsize)
{
	bin_size = binsize;
}

double CollimateParticleProcess::GetOutputBinSize() const
{
	return bin_size;
}

void CollimateParticleProcess::EnableImperfections(bool enable)
{
	Imperfections = enable;
}
} // end namespace ParticleTracking
