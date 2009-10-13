/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/18 17:52:43 $
// $Revision: 1.11 $
// 
/////////////////////////////////////////////////////////////////////////

#include "merlin_config.h"
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iterator>

#include "NumericalUtils/utils.h"
// CollimateParticleProcess
#include "BeamDynamics/ParticleTracking/CollimateParticleProcess.h"
// Aperture
#include "AcceleratorModel/Aperture.h"
// Spoiler
#include "AcceleratorModel/StdComponent/Spoiler.h"

using namespace std;

extern void ScatterParticle(PSvector& p, double X0, double x, double E0);

namespace {

using namespace ParticleTracking;

void OutputIndexParticles(const PSvectorArray lost_p, const list<size_t>& lost_i, ostream& os)
{
    PSvectorArray::const_iterator p = lost_p.begin();
    list<size_t>::const_iterator ip = lost_i.begin();

    while(p!=lost_p.end()) {
        os<<std::setw(12)<<right<<*ip;
        os<<*p;
        ++p;
        ++ip;
    }
}

}; // end anonymous namespace

namespace ParticleTracking {

CollimateParticleProcess::CollimateParticleProcess (int priority, int mode, std::ostream* osp)
        : ParticleBunchProcess("PARTICLE COLLIMATION",priority),cmode(mode),os(osp),
        createLossFiles(false),file_prefix(""),nstart(0),pindex(0),lossThreshold(1),scatter(false)
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

    if(active) {
        currentComponent = &component;
        s=0;
        Spoiler* aSpoiler = dynamic_cast<Spoiler*>(&component);
        is_spoiler = scatter && aSpoiler;

        if(!is_spoiler) { // not a spoiler so set up for normal hard-edge collimation
            at_entr = (COLL_AT_ENTRANCE & cmode)!=0;
            at_cent = (COLL_AT_CENTER & cmode)!=0;
            at_exit = (COLL_AT_EXIT & cmode)!=0;
            SetNextS();
        }
        else {
            at_entr=at_cent=false; // currently scatter only at exit
            at_exit = true;
            SetNextS();
            Xr = aSpoiler->GetMaterialRadiationLength();
			len = aSpoiler->GetLength();
        }
    }
    else {
        s_total += component.GetLength();
        currentComponent = 0;
    }
}

void CollimateParticleProcess::DoProcess (double ds)
{
    s+=ds;

    if(fequal(s,next_s)) {
        DoCollimation();
        SetNextS();
    }

    // If we are finished, GetNextS() will have set the process inactive.
    // In that case we can update s_total with the component length.
    if(!active)
        s_total += currentComponent->GetLength();
}

double CollimateParticleProcess::GetMaxAllowedStepSize () const
{
    return next_s-s;
}

void CollimateParticleProcess::IndexParticles (bool index)
{
    if(index && pindex==0)
        pindex = new list<size_t>;
    else if(!index && pindex!=0) {
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
    const Aperture *ap = currentComponent->GetAperture();

    PSvectorArray lost;
    list<size_t>  lost_i;

    list<size_t>::iterator ip;
    if(pindex!=0)
        ip=pindex->begin();

    for(PSvectorArray::iterator p = currentBunch->begin(); p!=currentBunch->end();) {
		if(!ap->PointInside((*p).x(),(*p).y(),s)) {

			// If the 'aperture' is a spoiler, then the particle is lost
			// if the DoScatter(*p) returns true (energy cut)
			if(!is_spoiler || DoScatter(*p)) {
				lost.push_back(*p);
				p=currentBunch->erase(p);
				if(pindex!=0) {
					lost_i.push_back(*ip);
					ip = pindex->erase(ip);
				}
			}
			else { // need to increment iterators
				p++;
				if(pindex!=0) {
					ip++;
				}
			}
		}
		else {
			p++;
			if(pindex!=0) {
				ip++;
			}
		}
	}

    nlost+=lost.size();
    DoOutput(lost,lost_i);

    if(double(nlost)/double(nstart)>=lossThreshold)
        throw ExcessiveParticleLoss(currentComponent->GetQualifiedName(),lossThreshold,nlost,nstart);
}

void CollimateParticleProcess::SetNextS ()
{
    if(at_entr) {
        next_s=0;
        at_entr=false;
    }
    else if(at_cent) {
        next_s=currentComponent->GetLength()/2;
        at_cent=false;
    }
    else if(at_exit) {
        next_s=currentComponent->GetLength();
        at_exit = false;
    }
    else
        active=false;
}

void CollimateParticleProcess::DoOutput (const PSvectorArray& lostb, const list<size_t>& lost_i)
{
    // Create a file and dump the lost particles
    // (if there are any)
    if(!lostb.empty()) {
        if(os!=0) {
            (*os)<<std::setw(24)<<left<<(*currentComponent).GetQualifiedName().c_str();
            (*os)<<std::setw(12)<<right<<currentComponent->GetLength();
            (*os)<<std::setw(12)<<right<<currentBunch->GetReferenceTime();
            (*os)<<std::setw(8)<<right<<lostb.size()<<endl;
        }
        if(createLossFiles) {
            string id = (*currentComponent).GetQualifiedName();
            pair<IDTBL::iterator,bool> result = idtbl.insert(IDTBL::value_type(id,0));
            int n = ++(*(result.first)).second;
            ostringstream fname;
            if(!file_prefix.empty())
                fname<<file_prefix;
            fname<<id<<'.'<<n<<".loss";
            ofstream file(fname.str().c_str());
            if(pindex==0)
                copy(lostb.begin(),lostb.end(), ostream_iterator<PSvector>(file) );
            else
                OutputIndexParticles(lostb,lost_i,file);
        }
    }
}

bool CollimateParticleProcess::DoScatter (Particle& p)
{
    double E0=currentBunch->GetReferenceMomentum();
    ScatterParticle(p,Xr,len,E0);
    return p.dp()<=-0.99; // return true if E below 1% cut-off
}

ExcessiveParticleLoss::ExcessiveParticleLoss (const string& c_id, double threshold, size_t nlost, size_t nstart)
        : MerlinException()
{
    ostringstream buffer;
    buffer<<"CollimateParticleProcess Exception\n";
    buffer<<"particle loss threshold of "<<100*threshold<<"% exceeded ";
    buffer<<'('<<nlost<<'/'<<nstart<<") at "<<c_id;
    SetMsg(buffer.str());
}

}; // end namespace ParticleTracking

