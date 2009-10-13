/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/10/24 19:15:24 $
// $Revision: 1.10 $
// 
/////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iterator>
#include <cassert>
#include "stdext/algorithm.h"

// ComponentFrame
#include "AcceleratorModel/Frames/ComponentFrame.h"
// TComponentFrame
#include "AcceleratorModel/Frames/TComponentFrame.h"
// Channels
#include "Channels/Channels.h"
// ring_iterator
#include "stdext/ring_iterator.h"
// ChannelServer
#include "AcceleratorModel/Implementation/ChannelServer.h"
// deleters
#include "stdext/deleters.h"
// AcceleratorModel
#include "AcceleratorModel/AcceleratorModel.h"
// SupportStructure
#include "AcceleratorModel/Supports/SupportStructure.h"

using namespace std;

namespace {

struct MatchName {
    StringPattern pattern;
    MatchName(const string& pat) : pattern(pat) {}
    bool operator()(const ComponentFrame* frm) {
        return pattern((*frm).GetComponent().GetQualifiedName());
    }
};

struct ModelStats {
    map<string,int>& s;
    ModelStats(map<string,int>& stats) :s(stats) {}
    void operator()(const ModelElement* element) {
        s[element->GetType()]++;
    }
};

class ExtractAcceleratorSupports : public FrameTraverser {
public:
	explicit ExtractAcceleratorSupports(AcceleratorSupportList& asList) 
		: asl(asList), nFound(0) {}
	void ActOn(LatticeFrame* frame);
	size_t nFound;
private:
	AcceleratorSupportList& asl;
};

void ExtractAcceleratorSupports::ActOn(LatticeFrame* aFrame)
{
	SupportStructure *ss = dynamic_cast<SupportStructure*>(aFrame);
	if(ss)
		nFound += ss->ExportSupports(asl);
}

} // end anonymous namespace

extern ChannelServer* ConstructChannelServer();

AcceleratorModel::AcceleratorModel ()
        : globalFrame(0)
{
    theElements = new ElementRepository();
    chServer = ConstructChannelServer();
    chServer->SetRepository(theElements);
}

AcceleratorModel::~AcceleratorModel ()
{
    if(chServer)
        delete chServer;
    if(theElements)
        delete theElements;
}

AcceleratorModel::Beamline AcceleratorModel::GetBeamline ()
{
    BeamlineIterator i = lattice.end();
    advance(i,-1);
    return Beamline(lattice.begin(),i,0,lattice.size()-1);
}

AcceleratorModel::Beamline AcceleratorModel::GetBeamline (AcceleratorModel::Index n1, AcceleratorModel::Index n2) throw (BadRange)
{
    if(n1<0 || n2>=lattice.size())
        throw BadRange();

    BeamlineIterator i1 = lattice.begin();
    BeamlineIterator i2 = lattice.begin();
    advance(i1,n1);
    advance(i2,n2);
    return Beamline(i1,i2,n1,n2);
}

AcceleratorModel::Beamline AcceleratorModel::GetBeamline (const string& pat1, const string& pat2, int n1, int n2) throw (BadRange)
{
    assert(n1>=1 && n2>=1);

    StringPattern p1(pat1),p2(pat2);
    BeamlineIterator i1=lattice.end();
    BeamlineIterator i2=lattice.end();
    int nn1(0),nn2(0);
    int ni=0, ni1, ni2;

    for(BeamlineIterator i = lattice.begin(); i!=lattice.end() && (nn1!=n1 || nn2!=n2); i++,ni++) {
        string id = (*i)->IsComponent() ? ((*i)->GetComponent()).GetQualifiedName() : (*i)->GetQualifiedName();
        if(nn1<n1 && p1(id) && ++nn1 == n1) {
            i1=i;
            ni1=ni;
        }
        if(nn2<n2 && p2(id) && ++nn2 == n2) {
            i2=i;
            ni2=ni;
        }
    }

    if(i1==lattice.end() || i2==lattice.end())
        throw BadRange();

    return Beamline(i1,i2,ni1,ni2);
}

AcceleratorModel::RingIterator AcceleratorModel::GetRing (int n)
{
    BeamlineIterator i = lattice.begin();
    advance(i,n);
    return RingIterator(lattice,i);
}

AcceleratorModel::Beamline AcceleratorModel::GetReversedBeamline ()
{
    BeamlineIterator i = lattice.end();
    advance(i,-1);
    return Beamline(i,lattice.begin(),lattice.size()-1,0);
}

int AcceleratorModel::ExtractComponents (const string& pat, vector<ComponentFrame*>& frames)
{
    vector<ComponentFrame*> results;
    if(pat=="*") // copy everything!
        copy(lattice.begin(),lattice.end(),back_inserter(results));
    else {
        MatchName mname(pat);
        for(BeamlineIterator bi=lattice.begin(); bi!=lattice.end(); bi++)
            if(mname(*bi))
                results.push_back(*bi);
    }
    frames.swap(results);
    return frames.size();
}

int AcceleratorModel::ExtractModelElements (const string& pat, vector<ModelElement*>& results)
{
    return theElements->Find(pat,results);
}

size_t AcceleratorModel::GetROChannels (const string& chID, std::vector<ROChannel*>& channels)
{
    return chServer->GetROChannels(chID,channels);
}

size_t AcceleratorModel::GetRWChannels (const string& chID, std::vector<RWChannel*>& channels)
{
    return chServer->GetRWChannels(chID,channels);
}

size_t AcceleratorModel::GetROChannels (AcceleratorModel::Beamline& aBeamline, const std::string& chid, std::vector<ROChannel*>& channels)
{
    return chServer->GetROChannels(aBeamline,chid,channels);
}

size_t AcceleratorModel::GetRWChannels (AcceleratorModel::Beamline& aBeamline, const std::string& chid, std::vector<RWChannel*>& channels)
{
    return chServer->GetRWChannels(aBeamline,chid,channels);
}

void AcceleratorModel::AddModelElement (ModelElement* element)
{
    if(element!=0)
        theElements->Add(element);
}

void AcceleratorModel::ReportModelStatistics (std::ostream& os) const
{
    using std::map;

    os<<"Arc length of beamline:     "<<globalFrame->GetGeometryLength()<<" meter"<<endl;
    os<<"Total number of components: "<<lattice.size()<<endl;
    os<<"Total number of elements:   "<<theElements->Size()<<endl;
    os<<endl;
    os<<"Model Element statistics\n";
    os<<"------------------------\n\n";

    map<string,int> stats;
    for_each(theElements->begin(),theElements->end(),ModelStats(stats));
    for(map<string,int>::iterator si=stats.begin(); si!=stats.end(); si++) {
        string atype = (*si).first;
        int count = (*si).second;
        os<<std::setw(20)<<left<<atype.c_str();
        os<<std::setw(4)<<count<<endl;
    }
    os<<endl;
}

size_t AcceleratorModel::GetIndecies(const std::string& pat,
                                     std::vector<AcceleratorModel::Index>& iarray) const
{
    return GetIndecies(const_cast<AcceleratorModel*>(this)->GetBeamline(),pat,iarray);
}

size_t AcceleratorModel::GetIndecies(const AcceleratorModel::Beamline& bline,
                                     const std::string& pat,
                                     std::vector<AcceleratorModel::Index>& iarray) const
{
    vector<Index> iarray1;
    StringPattern pattern(pat);
    
    Index n0 = distance(lattice.begin(),bline.begin());
    for(ConstBeamlineIterator fi = bline.begin(); fi!=bline.end(); fi++,n0++) {
        if((*fi)->IsComponent()) {
            string id = (*(*fi)).GetComponent().GetQualifiedName();
            if(pattern(id))
                iarray1.push_back(n0);
        }
    }
    iarray.swap(iarray1);
    return iarray.size();
}

size_t AcceleratorModel::GetAcceleratorSupports(AcceleratorSupportList& supports)
{
	ExtractAcceleratorSupports eas(supports);
	globalFrame->Traverse(eas);
	return eas.nFound;
}
