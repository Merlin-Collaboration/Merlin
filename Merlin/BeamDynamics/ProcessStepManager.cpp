/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iomanip>
#include "NumericalUtils/utils.h"
#include "BeamDynamics/BunchProcess.h"
#include "BeamDynamics/ProcessStepManager.h"
#include "AcceleratorModel/AcceleratorComponent.h"
#include "stdext/deleters.h"

namespace {

#ifndef  _MSC_VER
#define _MIN(a,b) std::min(a,b)
#define _MAX(a,b) std::max(a,b)
#endif

using namespace std;
using std::setw;

typedef list<BunchProcess*>::iterator proc_itor;
typedef list<BunchProcess*>::const_iterator const_proc_itor;
typedef list<BunchProcess*>::reverse_iterator rev_proc_itor;
typedef list<BunchProcess*>::const_reverse_iterator const_rev_proc_itor;

struct InitProc {
    Bunch& ibunch;
    InitProc(Bunch& b) : ibunch(b) {}
    void operator()(BunchProcess* proc) {
    	//cout<<"Debug: "<<proc->GetID()<<end;;
    	proc->InitialiseProcess(ibunch);
    }
};

struct SetCmpnt {
    AcceleratorComponent& cmp;
    SetCmpnt(AcceleratorComponent& ac) : cmp(ac) {}
    void operator()(BunchProcess* proc) {
        proc->SetCurrentComponent(cmp);
    }
};

struct CalcStepSize {
    double ds;
    CalcStepSize(double s_max) : ds(s_max) {}
    void operator()(BunchProcess* proc) {
        if(proc->IsActive())
            ds = _MIN(ds,proc->GetMaxAllowedStepSize());
    }
};

struct DoProc {
    double s0;
    double ds;
    const string& cid;
    ostream* vos;

    DoProc(double s, double ds1, const string& id, ostream* os)
            : s0(s), ds(ds1), cid(id), vos(os) {}

    void operator()(BunchProcess* proc) {
        if(proc->IsActive()) {
            if(vos!=0) Trace(proc);
            proc->DoProcess(ds);
        }
    }

    void Trace(BunchProcess*);
};

void DoProc::Trace(BunchProcess* proc)
{
    (*vos)<<setw(24)<<left<<cid.c_str();
    (*vos)<<setw(24)<<left<<(*proc).GetID().c_str();
    (*vos)<<"from: "<<right<<s0<<" to: "<<s0+ds<<" (step = "<<ds<<")"<<endl;
}

}; // end of anonymous namespace


ProcessStepManager::ProcessStepManager ()
        : total_s(0),log(0),processTable()
{}

ProcessStepManager::~ProcessStepManager ()
{
    for_each(processTable.begin(),processTable.end(),deleter<BunchProcess>());
}

void ProcessStepManager::Initialise (Bunch& bunch)
{
    for_each(processTable.begin(),processTable.end(),InitProc(bunch));
    total_s=0;
}

void ProcessStepManager::Track (AcceleratorComponent& component)
{
    const string id = component.GetQualifiedName();

    for_each(processTable.begin(),processTable.end(),SetCmpnt(component));

    const double sc = component.GetLength();
    double s=0;
    do{
        double ds = for_each(processTable.begin(),processTable.end(),CalcStepSize(sc-s)).ds;
        for_each(processTable.begin(),processTable.end(),DoProc(s,ds,id,log));
        s+=ds;
    } while(!fequal(sc,s));

    total_s += sc;
}

double ProcessStepManager::GetIntegratedLength ()
{
    return total_s;
}

void ProcessStepManager::AddProcess (BunchProcess* aProcess)
{
    if(aProcess->GetPriority()<0)
        processTable.push_back(aProcess);
    else {
        proc_itor p;
        for(p = processTable.begin();
                p!=processTable.end() && (*p)->GetPriority() <= aProcess->GetPriority();
                p++);
        processTable.insert(p,aProcess);
    }
}

bool ProcessStepManager::RemoveProcess (BunchProcess* aProcess)
{
    proc_itor p = find(processTable.begin(),processTable.end(),aProcess);
    processTable.erase(p);
    return p!=processTable.end();
}

void ProcessStepManager::ClearProcesses ()
{
    for_each(processTable.begin(),processTable.end(),deleter<BunchProcess>());
    processTable.clear();
}

void ProcessStepManager::SetLogStream (ostream* os)
{
    log=os;
}

