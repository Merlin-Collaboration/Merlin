/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/20 13:42:54 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#include "merlin_config.h"
#include <cassert>
#include <map>
#include <iomanip>
#include <exception>
#include <typeinfo>

// Drift
#include "AcceleratorModel/StdComponent/Drift.h"
// SequenceFrame
#include "AcceleratorModel/Frames/SequenceFrame.h"
// AcceleratorModelConstructor
#include "AcceleratorModel/Construction/AcceleratorModelConstructor.h"

namespace {

struct ModelStats {

    map<string,int>& s;

    ModelStats(map<string,int>& stats) :s(stats) {}

    void operator()(const ModelElement* element)
    {
        s[element->GetType()]++;
    }

};

// class MEExtractor
// FrameTraverser implementation used by AppendFrame(SequenceFrame*)
// to add the specified frame to the accelerator model.

class MEExtractor : public FrameTraverser {
public:

    MEExtractor(ElementRepository* erp, AcceleratorModel::FlatLattice* lp)
            : erepo(erp),lattice(lp)
    {}

    virtual void ActOn(LatticeFrame* frame) {
        erepo->Add(frame);
        try {
            ComponentFrame* cf = dynamic_cast<ComponentFrame*>(frame);
            if(cf) {
                erepo->Add(&(cf->GetComponent()));
                lattice->push_back(cf);
            }
        }
        catch(bad_cast& b) {
            cerr<<"bad cast cought here!"<<endl;
        }
        catch(...) {
            cerr<<"something weird here"<<endl;
        }
    };

private:
    ElementRepository* erepo;
    AcceleratorModel::FlatLattice* lattice;
};


}; // end anonymous namespace


AcceleratorModelConstructor::AcceleratorModelConstructor ()
        : currentModel(0)
{}

AcceleratorModelConstructor::~AcceleratorModelConstructor ()
{
    if(currentModel)
        delete currentModel;
}

void AcceleratorModelConstructor::NewModel ()
{
    if(currentModel) {
        delete currentModel;
        while(!frameStack.empty())
            frameStack.pop();
    }

    currentModel = new AcceleratorModel();
    SequenceFrame *globalFrame = new SequenceFrame("GLOBAL",SequenceFrame::originAtEntrance);
    currentModel->globalFrame = globalFrame;
    frameStack.push(globalFrame);
}

AcceleratorModel* AcceleratorModelConstructor::GetModel ()
{
    frameStack.pop();
    assert(currentModel && frameStack.empty());

    currentModel->globalFrame->ConsolidateConstruction();

    AcceleratorModel* t=currentModel;
    currentModel=0;
    return t;
}

void AcceleratorModelConstructor::NewFrame (SequenceFrame* aFrame)
{
    assert(currentModel && aFrame);
    currentModel->theElements->Add(aFrame);
    frameStack.push(aFrame);
}

void AcceleratorModelConstructor::EndFrame ()
{
    SequenceFrame* cf = frameStack.top();
    frameStack.pop();
    (frameStack.top())->AppendFrame(*cf);
}

void AcceleratorModelConstructor::AppendDrift (double d)
{
    Drift* newd = new Drift("UNNAMED",d);
    AppendComponentFrame(new TComponentFrame<Drift>(*newd));
}

void AcceleratorModelConstructor::AddModelElement (ModelElement* element)
{
    currentModel->theElements->Add(element);
}

void AcceleratorModelConstructor::ReportStatistics (std::ostream& os) const
{
    using std::map;

    os<<"Arc length of beamline:     "<<currentModel->globalFrame->GetGeometryLength()<<" meter"<<endl;
    os<<"Total number of components: "<<(*currentModel).lattice.size()<<endl;
    os<<"Total numner of elements:   "<<currentModel->theElements->Size()<<endl;
    os<<endl;
    os<<"Model Element statistics\n";
    os<<"------------------------\n\n";

    map<string,int> stats;
    for_each(currentModel->theElements->begin(),currentModel->theElements->end(),ModelStats(stats));
    for(map<string,int>::iterator si=stats.begin(); si!=stats.end(); si++) {
        string atype = (*si).first;
        int count = (*si).second;
        os<<setw(20)<<left<<atype.c_str();
        os<<setw(4)<<count<<endl;
    }
    os<<endl;
}

void AcceleratorModelConstructor::AppendFrame(SequenceFrame* aFrame)
{
    assert(currentModel);

    // We extract all the model elements from aFrame
    // and add them to the repository. In addition, we
    // extract (in order) all the ComponentFrame objects
    // and append them to the lattice.

    MEExtractor mextr(currentModel->theElements,&(currentModel->lattice));
    aFrame->Traverse(mextr);

    (frameStack.top())->AppendFrame(*aFrame);
}

void AcceleratorModelConstructor::AppendComponentFrame (ComponentFrame* cf)
{
    assert(currentModel);
    currentModel->theElements->Add(cf);
    if(cf->IsComponent())
        currentModel->theElements->Add(&(cf->GetComponent()));
    (*currentModel).lattice.push_back(cf);
    (frameStack.top())->AppendFrame(*cf);
	cf->SetBeamlineIndex((*currentModel).lattice.size()-1);
}

