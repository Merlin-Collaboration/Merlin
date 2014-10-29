/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
// 
/////////////////////////////////////////////////////////////////////////

#include "merlin_config.h"
#include <iterator>
#include <iomanip>
#include <typeinfo>
#include <fstream>
#include <sstream>
#include "AcceleratorModel/Apertures/CollimatorAperture.h"
#include "NumericalUtils/utils.h"
#include "NumericalUtils/PhysicalUnits.h"
// CollimateParticleProcess
#include "Collimators/CollimateParticleProcess.h"
// Aperture
#include "AcceleratorModel/Aperture.h"
// Spoiler
#include "AcceleratorModel/StdComponent/Spoiler.h"

#include "BeamDynamics/ParticleTracking/ParticleComponentTracker.h"

#include "AcceleratorModel/Apertures/InterpolatedApertures.h"

using namespace std;

//extern void ScatterParticle(PSvector& p, double X0, double x, double E0);
//extern void ScatterProton(PSvector& p, double x, double E0, const TiltedAperture* tap);

namespace {

using namespace ParticleTracking;

void OutputIndexParticles(const PSvectorArray lost_p, const list<size_t>& lost_i, ostream& os)
{
    PSvectorArray::const_iterator p = lost_p.begin();
    list<size_t>::const_iterator ip = lost_i.begin();

    while(p!=lost_p.end())
    {
        os<<std::setw(12)<<right<<*ip;
        os<<*p;
        ++p;
        ++ip;
    }
}

} // end anonymous namespace

namespace ParticleTracking {

CollimateParticleProcess::CollimateParticleProcess (int priority, int mode, std::ostream* osp)
        : ParticleBunchProcess("PARTICLE COLLIMATION",priority),cmode(mode),os(osp),
        createLossFiles(false),file_prefix(""),lossThreshold(1),nstart(0),pindex(0),scatter(false),bin_size(0.1*PhysicalUnits::meter),Imperfections(false)
{}

CollimateParticleProcess::~CollimateParticleProcess ()
{
    if(pindex!=0)
        delete pindex;
}

void CollimateParticleProcess::InitialiseProcess (Bunch& bunch)
{
    ParticleBunchProcess::InitialiseProcess(bunch);
    idtbl.clear();
    if(currentBunch) {
        nstart = currentBunch->size();
        nlost = 0;
        if(pindex!=0) {
            pindex->clear();
            for(size_t n=0; n<nstart; n++)
                pindex->push_back(n);
        }
    }
}

void CollimateParticleProcess::SetCurrentComponent (AcceleratorComponent& component)
{
	active = (currentBunch!=0) && (component.GetAperture()!=0);
	if(active)
	{
//		cout << endl;
		currentComponent = &component;
//		cout << currentComponent->GetQualifiedName() << "\t" << nstart << endl;
//		cout << currentComponent->GetQualifiedName() << "\t" << currentComponent->GetAperture()->GetApertureType() << endl;
		s=0;
		Spoiler* aSpoiler = dynamic_cast<Spoiler*>(&component);
		is_spoiler = scatter && aSpoiler;

		if(!is_spoiler)
		{ // not a spoiler so set up for normal hard-edge collimation
			at_entr = (COLL_AT_ENTRANCE & cmode)!=0;
			at_cent = (COLL_AT_CENTER & cmode)!=0;
			at_exit = (COLL_AT_EXIT & cmode)!=0;
			SetNextS();
		}
		else
		{
			//at_entr = at_cent = false; // currently scatter only at exit
			at_cent = false; // currently scatter only at exit
			at_exit = at_entr = true;
			SetNextS();
			//Xr = aSpoiler->GetMaterialRadiationLength();
			currentBunch->SetScatterConfigured(false);
			len = aSpoiler->GetLength();
		//	CollimatorAperture* CollimatorJaw = dynamic_cast<CollimatorAperture*>(aSpoiler->GetAperture());
		}

		//For precision tracking of lost particles in non-collimators
		//make a copy of the input array if we are a magnet.
		//This should also occur before any tracking, hence also any poleface rotations on a dipole, etc
		if(!is_spoiler)
		{
			InputArray = currentBunch->GetParticles();
		}
	}
	else
	{
		s_total += component.GetLength();
		currentComponent = 0;
	}
}

void CollimateParticleProcess::DoProcess (double ds)
{
    //std::cout << "Call DoProcess " << std::endl;
    s+=ds;

    if(fequal(s,next_s))
    {
	//This lets the scattering routine know how far down the collimator we are for aperture checking inside the scattering step.
	currentBunch->SetIntS(s-ds);
        DoCollimation();
        SetNextS();
    }

    // If we are finished, GetNextS() will have set the process inactive.
    // In that case we can update s_total with the component length.
    if(!active)
        s_total += currentComponent->GetLength();
	//cout<<"the component name is:"<<currentComponent->GetName()<<endl;
	//cout<<"component name:"<<currentComponent->GetQualifiedName()<<endl;
	//cout << s << "\t" << s_total << "\t" << currentComponent->GetLength() << endl;
}

double CollimateParticleProcess::GetMaxAllowedStepSize () const
{
	if(!is_spoiler)
	{
		return next_s-s;
	}
	else
	{
		return bin_size;
	}
}

void CollimateParticleProcess::IndexParticles (bool index)
{
	if(index && pindex==0)
	{
		pindex = new list<size_t>;
	}
	else if(!index && pindex!=0)
	{
		delete pindex;
		pindex=0;
	}
}

void CollimateParticleProcess::IndexParticles (list<size_t>& anIndex)
{
    if(!pindex)
        delete pindex;

    pindex=&anIndex;
}


void CollimateParticleProcess::SetLossThreshold (double losspc)
{
	lossThreshold = losspc/100.0;
}

void CollimateParticleProcess::DoCollimation ()
{
	//std::cout << "Call DoCollimation " << std::endl;
//	cout << "CurrentBunch size: " << currentBunch->size() << endl;
	//The apture of this element
	const Aperture *ap = currentComponent->GetAperture();

	//The array of lost particles
	PSvectorArray lost;
	list<size_t>  lost_i;

	list<size_t>::iterator ip;
	if(pindex!=0)
	ip=pindex->begin();

	//For copying surviving particles to, which is faster than deleting the individual lost particles
	ParticleBunch* NewBunch=new ParticleBunch(currentBunch->GetReferenceMomentum(),currentBunch->GetTotalCharge()/currentBunch->size());

/*
	//For precision tracking of lost particles in non-collimators
	PSvectorArray InputArray;				//The input array
	std::vector<unsigned int> LostParticlePositions;	//A list of particles we want to use in the input array

	//make a copy of the input array if we are a magnet.
	if(!is_spoiler)
	{
		InputArray = currentBunch->GetParticles();
	}
*/

/*	cout << currentComponent->GetQualifiedName() << "\t" << currentComponent->GetLength() << "\t" << s << endl;
	if(currentComponent->GetQualifiedName() != "Spoiler.TCP.C6L7.B1")
	{
		abort();
	}
*/

	unsigned int particle_number=0;

	if(is_spoiler) 
	{
		for(PSvectorArray::iterator p = currentBunch->begin(); p!=currentBunch->end();)
		{
			(*p).x() -= bin_size * (*p).xp();
			(*p).y() -= bin_size * (*p).yp();
			p++;
		}
	}

	for(PSvectorArray::iterator p = currentBunch->begin(); p!=currentBunch->end();)
	{
		//If we are collimating at the end of the element, track back a drift
		//Do not do this at the start of the element.

//		if(is_spoiler) 
//		{
//			(*p).x() -= bin_size * (*p).xp();
//			(*p).y() -= bin_size * (*p).yp();
//		}

//		if(!ap->PointInside(( *p).x(), (*p).y(), s + (*p).ct() ))
		if(!ap->PointInside( (*p).x(), (*p).y(), s) )
		{
			// If the 'aperture' is a spoiler, then the particle is lost
			// if the DoScatter(*p) returns true (energy cut)
			// If not a spoiler, then do not scatter and directly remove the particle.
			if(!is_spoiler || DoScatter(*p))
			{
//				cout << "Lost Particle at: " << s + (*p).ct() << endl;
				if(is_spoiler)
				{
					(*p).ct() += (s-bin_size);
				}

				lost.push_back(*p);
				/* This is slow for a STL Vector - instead we place the surviving particles into a new bunch and then swap - this is faster */
				//p=currentBunch->erase(p);
				p++;
				if(pindex!=0)
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
				if(pindex!=0)
				{
					ip++;
				}
			}
		}
		else
		{
			if(is_spoiler) 
			{
				(*p).x() += bin_size * (*p).xp();
				(*p).y() += bin_size * (*p).yp();
			}

			//Not interacting with the collimator: "Inside" the aperture; particle lives
			NewBunch->AddParticle(*p);
			p++;
			if(pindex!=0)
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
if(LostParticlePositions.size() != 0 && !is_spoiler)
{
//cout << currentComponent->GetQualifiedName() << "\tLostParticlePositions: " << LostParticlePositions.size() << endl;
	//make a new particle bunch to track the lost particles
	ParticleBunch* LostBunch=new ParticleBunch(currentBunch->GetReferenceMomentum(),currentBunch->GetTotalCharge()/currentBunch->size());
	double length = currentComponent->GetLength();
	//If we are dealing with a non-zero length element, we must do tracking
	if(length != 0)
	{
		//Clear out the old lost particles, these will all be at the end of the element which we do not want for a magnet.
		lost.clear();

		//Grab the lost particles from the copied input particle array and add them to the new particle bunch
		for(vector<unsigned int>::iterator p = LostParticlePositions.begin(); p!=LostParticlePositions.end(); p++)
//		for(int i =0; i < InputArray.size(); i++)
		{
			LostBunch->AddParticle(InputArray[*p]);
//			LostBunch->AddParticle(InputArray[i]);
		}
//		cout << "LostBunch size: " << LostBunch->size() << "\tLost size: " << lost.size() << endl;
		//Create a new tracker
		ParticleComponentTracker* LostParticleTracker = new ParticleComponentTracker();

		//Tell the new tracker to use the particle bunch of lost particles that has been made
		LostParticleTracker->SetBunch(*LostBunch);

		//Prepare the tracker to use the current accelerator component
		currentComponent->PrepareTracker(*LostParticleTracker);
//		cout << LostBunch->size() << endl;
		//While we are still inside the component
		while((LostParticleTracker->GetRemainingLength() ) >= 0 && LostBunch->size() != 0)
		{
			double StepSize = bin_size;
			//If the remaining length of component is less than the step size, set the step size to this value
			if((LostParticleTracker->GetRemainingLength() - StepSize) < 0)
			{
				StepSize = LostParticleTracker->GetRemainingLength();
			}

			//Track the appropriate step length
			/*cout << "Tracking: " << currentComponent->GetQualifiedName() << "\tStepsize: " << StepSize << "\tPreStep: " << LostParticleTracker->GetIntegratedLength() << \
				"\tPostStep: " << LostParticleTracker->GetIntegratedLength() + StepSize << "\tLength: " << currentComponent->GetLength() << endl;*/

			double IntegratedLength = LostParticleTracker->GetIntegratedLength();
			//Now loop over each particle in turn
			for(PSvectorArray::iterator p = LostBunch->begin(); p!=LostBunch->end();)
			{
				//Check if the particle is outside the aperture
				//s, is where the integrator will start
				//LostParticleTracker->GetIntegratedLength() will give the position integrated past this point
				//(*p).ct() will give the offset for this specific particle

//				if(!ap->PointInside((*p).x(),(*p).y(),IntegratedLength + (*p).ct() ))
				if(!ap->PointInside((*p).x(),(*p).y(),IntegratedLength ))
				{
//					cout << "Lost Particle at: " << IntegratedLength  << endl;
					//if not, delete the particle, and add the coordintes to the lost bunch list (PSvectorArray lost)
					/*cout <<	"Lost at: " << LostBunch->size() << "\t" << currentComponent->GetQualifiedName() << "\t" << \
					LostParticleTracker->GetIntegratedLength() << "\t" << length << endl;*/

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
					p=LostBunch->erase(p);
				}
				//else, the particle is inside and can be kept for this step
				else
				{
					p++;
				}
			}

			//Now move forward...
			if((LostParticleTracker->GetRemainingLength() ) > 0)
			{
				LostParticleTracker->TrackStep(StepSize);			
			}
		}
//		cout << "LostBunch size: " << LostBunch->size() << "\tLost size: " << lost.size() << endl;

		//If there is anything left - possible bug.
		if(LostBunch->size() != 0)
		{
//			abort();
			cout <<	"POSSIBLE BUG: Leftovers: " << LostBunch->size() << "\t" << currentComponent->GetQualifiedName() << "\t" << \
			LostParticleTracker->GetIntegratedLength() << "\t" << length << endl;
			for(PSvectorArray::iterator p = LostBunch->begin(); p!=LostBunch->end(); p++)
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
	//if the element has zero length nothing needs to be done since all the losses will have occured at the same point anyway.
	//So PSvectorArray loss will contain the correct information
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//cout << "Final lost size:" << lost.size() << endl;
//cout << "Final CurrentBunch size: " << currentBunch->size() << endl;
	nlost+=lost.size();
	//cout << currentComponent->GetQualifiedName() << "\t" << nlost << "\t" << nstart << "\t" << currentBunch->size() << endl;
	//cout << "The number of particles lost is: " << nlost << endl;
	DoOutput(lost,lost_i);

	//make sure to clear up
	InputArray.clear();				//The input array
	LostParticlePositions.clear();	//A list of particles we want to use in the input array

	if(double(nlost)/double(nstart)>=lossThreshold)
	{
		cout << "nlost: " << nlost << "\tnstart: " << nstart << endl;
		throw ExcessiveParticleLoss(currentComponent->GetQualifiedName(),lossThreshold,nlost,nstart);
	}
}


void CollimateParticleProcess::SetNextS ()
{

	if(at_entr)
	{
		next_s = 0;
		at_entr = false;
	}
	else if(at_cent)
	{
		next_s=currentComponent->GetLength()/2;
		at_cent=false;
	}
	else if(at_exit)
	{
		next_s=currentComponent->GetLength();
		at_exit = false;
	}
	else
	{
		active=false;
	}

	if(is_spoiler)
	{
		active = true;
		next_s = s + bin_size;
	}
}

void CollimateParticleProcess::DoOutput (const PSvectorArray& lostb, const list<size_t>& lost_i)
{
        //cout << "Call Do output" << endl;
	// Create a file and dump the lost particles
	// (if there are any)
//	cout << currentComponent->GetQualifiedName() << endl;
	if(!lostb.empty())
	{
//	PSvectorArray lostp = bin_lost_output(lostb);
        	if(os!=0)
		{
			double length = currentComponent->GetLength();
			double** lostp;
			int n;
			//deal with zero length elements separately
			if(length != 0)
			{
//			cout << "non zero length" << endl;

				//We bin in 0.1m sections of the element, this is how many bins we have
				//n is our number of bins
				n = (length / bin_size) +1;

				bool overflow = false;

				//The element length may not be an exact multiple of the bin size, so we must take this into account.
				if ((length/(double)n) != 0.0)
				{
					//Add an aditional bin
					n++;
					overflow = true;
				}
//			cout << n << "\t" << length << "\t" << bin_size << "\t" << overflow << endl;
				//Create the array - n rows, 3 cols
				lostp = new double*[n];
				for(int i = 0; i<n; ++i)
				{
					lostp[i] = new double[3];
				}

				//Initialize the array elements
				for(int j = 0; j<n; j++)
				{
					lostp[j][0] = bin_size*j;	//Start position
					lostp[j][1] = bin_size;		//Length
					lostp[j][2] = 0.0;		//Entries
				}

				//deal with the last bin if needed - will have a different length.
				if(overflow)
				{
					lostp[n-1][0] = bin_size*(n-1);
					lostp[n-1][1] = length - (bin_size*(n-1));
					lostp[n-1][2] = 0.0;
				}

				//cout << "Number of lost particles at " << currentComponent->GetName() << ": " << lostb.size() << endl;
//				cout << "I bet it crashes here: " << lostb.size() << endl;
				for(size_t l=0; l<lostb.size(); l++)
				{
					int x = (lostb[l].ct())/bin_size;
//					cout << lostb[l].ct() << "\t" << bin_size << "\t" << x << endl;
					//Add one loss count to this bin
					lostp[x][2]++;
				}
			}
			//now deal with zero length elements;
			else
			{
//			cout << "zero length" << endl;
				n = 1;
				lostp = new double*[1];
				lostp[0] = new double[3];
				lostp[0][0] = 0.0;		//Start position
				lostp[0][1] = bin_size;		//Length
				lostp[0][2] = lostb.size();	//Entries
			}

			//cout << "output bin loop" << endl;
			//Now to do the output - first loop over each bin
			for(int j = 0; j<n; j++)
			{
				//We then check if there are any lost particles in this bin
				//If there are lost particles, then we must output the count (lost count > 0)
				if(lostp[j][2] > 0.0)
				{
					//cout << lostp[j][0] << "\t" << lostp[j][1] << "\t" << lostp[j][2] << endl;
					(*os).precision(16);
					(*os) << std::setw(35)<<left<<(*currentComponent).GetQualifiedName().c_str();			//Component name - can be quite long for certain LHC magnets
					//(*os) << std::setw(24)<<left<<lostp[j][0] + currentBunch->GetReferenceTime()-length;	//Bin start position
					(*os) << std::setw(24)<<left<<lostp[j][0] + currentComponent->GetComponentLatticePosition();	//Bin start position
					(*os) << std::setw(24)<<left<<lostp[j][1];							//Bin length
					(*os) << "\t" << lostp[j][2];									//Loss count
					(*os) << "\t" << j * lostp[0][1];								//Loss position in element
					(*os) << endl;
				}
			}
			delete [] lostp;
		/*
			(*os)<<std::setw(24)<<left<<(*currentComponent).GetQualifiedName().c_str();
			(*os)<<std::setw(12)<<right<<currentComponent->GetLength();
			(*os)<<std::setw(12)<<right<<currentBunch->GetReferenceTime();
			(*os)<<std::setw(8)<<right<<lostb.size()<<endl;
		*/
		
		
		}
		//cout << "createlossfiles" << endl;
		if(createLossFiles)
		{
			string id = (*currentComponent).GetQualifiedName();
			pair<IDTBL::iterator,bool> result = idtbl.insert(IDTBL::value_type(id,0));
			int n = ++(*(result.first)).second;
			ostringstream fname;

			if(!file_prefix.empty())
			{
				fname << file_prefix;
			}

			fname << "_" << currentComponent->GetComponentLatticePosition() << "_";
			fname << id << '.' << n << ".loss";
			ofstream file(fname.str().c_str());
			if (!file) {cerr << "CollimateParticleProcess::DoOutput(): Failed to open " << fname.str() << endl; exit(1);}

			if(pindex==0)
			{
				copy(lostb.begin(),lostb.end(),ostream_iterator<PSvector>(file));
			}
			else
			{
				OutputIndexParticles(lostb,lost_i,file);
			}
		}
	}
//cout << "end Do output" << endl;
}

bool CollimateParticleProcess::DoScatter (Particle& p)
{
	//std::cout << "DoScatter " << std::endl;
	//double E0=currentBunch->GetReferenceMomentum();
	const CollimatorAperture *tap = (CollimatorAperture*) currentComponent->GetAperture();

	//int scatter_type = currentBunch->Scatter(p,len,tap);
	int scatter_type = currentBunch->Scatter(p,bin_size,tap);

//	return p.dp()<=-0.99; // return true if E below 1% cut-off
	if(scatter_type == 1)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

ExcessiveParticleLoss::ExcessiveParticleLoss (const string& c_id, double threshold, size_t nlost, size_t nstart)
        : MerlinException()
{
	ostringstream buffer;
	buffer << "CollimateParticleProcess Exception\n";
	buffer << "particle loss threshold of " << 100*threshold << "% exceeded ";
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



//cout << "LostBunch:\t" << LostBunch->size() << endl;
//					cout << "not lost at: " << s+(*p).ct() << endl;
//				cout << (*p).x() << "\t" << (*p).y() << "\t" << s << "\t" << LostParticleTracker->GetRemainingLength() << "\t" << LostBunch->size() << endl;
//				cout << "delete: " << LostBunch->size() << endl;
//				cout << "lost at: " << s+(*p).ct() << endl;
//cout << "END: " << currentComponent->GetQualifiedName() << endl;

//cout << endl <<"START: " << currentComponent->GetQualifiedName() << endl;
//cout << "Input size:\t" << InputArray.size() << endl;
//cout << "Lost size:\t" << LostParticlePositions.size() << endl;
//cout << "Lost!: " << currentComponent->GetQualifiedName() << "\t" << currentComponent->GetLength() << "\t" << currentComponent->GetComponentLatticePosition() << endl;
//cout << "Lost!: " << currentComponent->GetQualifiedName() << "\t" << foo->GetIntegratedLength() << endl;
//cout << currentComponent->GetQualifiedName() << endl;
//cout << "LOSTBUNCH SIZE: " << LostBunch->size() << endl;
//	cout << "New Tracker" << endl;
//ParticleBunchProcess::InitialiseProcess(bunch);
//	cout << "prep Tracker" << endl;
//	cout << LostParticleTracker->GetState() << endl;
//	cout << "check loss " << LostParticleTracker->GetRemainingLength() << endl;
//		cout << "Using step size: " << StepSize << endl;
//		cout << "TrackStep " << endl;
//cout << "LOSTBUNCH SIZE: " << LostBunch->size() << endl;

/*void CollimateParticleProcess::bin_lost_output(const PSvectorArray& lostb)
{

	double leng = currentComponent->GetLength();
	int n = leng / 0.1;


	bool overflow = false;
	if ((leng/(double)n) != 0.0)
	{
		n++;
		overflow = true;
	}


	//Create the array
	double** lostp = new double*[n];
	for(int i = 0; i<n; ++i)
	{
		lostp[i] = new double[3];
	}
	//fill the aray
	for(int j = 0; j<n; j++)
	{
		lostp[j][0] = 0.1*j;	//Start position
		lostp[j][1] = 0.1;	//Length
		lostp[j][2] = 0.0;	//Entries
	}
	//deal with the last bin
	if(overflow)
	{
		lostp[n-1][0] = 0.1*(n-1);
		lostp[n-1][1] = leng - (0.1*(n-1));
		lostp[n-1][2] = 0.0;
	}


	for(size_t l = 0; l < lostb.size(); l++)
	{
		int x = lostb[l].ct()/0.1;
		lostp[x][2]++;
	}

	for(int j = 0; j<n; j++)
	{
		cout << lostp[j][0] << "\t" << lostp[j][1] << "\t" << lostp[j][2] << endl;
	}
}
*/
