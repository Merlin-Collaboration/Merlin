/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <algorithm>
#include <iomanip>
#include "utils.h"
#include "BunchProcess.h"
#include "ProcessStepManager.h"
#include "AcceleratorComponent.h"
#include "deleters.h"

#include "MerlinProfile.h"

namespace
{

#ifndef  _MSC_VER
#define _MIN(a, b) std::min(a, b)
#define _MAX(a, b) std::max(a, b)
#endif

using std::setw;

typedef std::list<BunchProcess*>::iterator proc_itor;
typedef std::list<BunchProcess*>::const_iterator const_proc_itor;
typedef std::list<BunchProcess*>::reverse_iterator rev_proc_itor;
typedef std::list<BunchProcess*>::const_reverse_iterator const_rev_proc_itor;

struct InitProc
{
	Bunch& ibunch;
	InitProc(Bunch& b) :
		ibunch(b)
	{
	}
	void operator()(BunchProcess* proc)
	{
		MERLIN_PROFILE_START_TIMER(proc->GetID());
		proc->InitialiseProcess(ibunch);
		MERLIN_PROFILE_END_TIMER(proc->GetID());
	}

};

struct SetCmpnt
{
	AcceleratorComponent& cmp;
	SetCmpnt(AcceleratorComponent& ac) :
		cmp(ac)
	{
	}
	void operator()(BunchProcess* proc)
	{
		proc->SetCurrentComponent(cmp);
	}

};

struct CalcStepSize
{
	double ds;
	CalcStepSize(double s_max) :
		ds(s_max)
	{
	}
	void operator()(BunchProcess* proc)
	{
		if(proc->IsActive())
		{
			ds = _MIN(ds, proc->GetMaxAllowedStepSize());
		}
	}

};

struct DoProc
{
	double s0;
	double ds;
	const string& cid;
	ostream* vos;

	DoProc(double s, double ds1, const string& id, ostream* os) :
		s0(s), ds(ds1), cid(id), vos(os)
	{
	}

	void operator()(BunchProcess* proc)
	{
		if(proc->IsActive())
		{
			if(vos != nullptr)
			{
				Trace(proc);
			}
			MERLIN_PROFILE_START_TIMER(proc->GetID());
//	    cout << proc->GetID() << "\t" << ds << endl;
			proc->DoProcess(ds);
			MERLIN_PROFILE_END_TIMER(proc->GetID());
		}
	}

	void Trace(BunchProcess*);

};

void DoProc::Trace(BunchProcess* proc)
{
	(*vos) << setw(30) << left << cid.c_str();
	(*vos) << setw(24) << left << (*proc).GetID().c_str();
	(*vos) << "from: " << right << s0 << " to: " << s0 + ds << " (step = " << ds << ")" << std::endl;
}

} // end of anonymous namespace

ProcessStepManager::ProcessStepManager() :
	total_s(0), log(nullptr), processTable()
{
}

ProcessStepManager::~ProcessStepManager()
{
	for_each(processTable.begin(), processTable.end(), deleter<BunchProcess>());
}

void ProcessStepManager::Initialise(Bunch& bunch)
{
	for_each(processTable.begin(), processTable.end(), InitProc(bunch));
	total_s = 0;
}

void ProcessStepManager::Track(AcceleratorComponent& component)
{
	const std::string id = component.GetQualifiedName();

	for_each(processTable.begin(), processTable.end(), SetCmpnt(component));

	const double sc = component.GetLength();
	double s = 0;
	do
	{
		double ds = for_each(processTable.begin(), processTable.end(), CalcStepSize(sc - s)).ds;
		for_each(processTable.begin(), processTable.end(), DoProc(s, ds, id, log));
		s += ds;
	} while(!fequal(sc, s));

	total_s += sc;
}

double ProcessStepManager::GetIntegratedLength()
{
	return total_s;
}

void ProcessStepManager::AddProcess(BunchProcess* aProcess)
{
	if(aProcess->GetPriority() < 0)
	{
		processTable.push_back(aProcess);
	}
	else
	{
		proc_itor p;
		for(p = processTable.begin();
			p != processTable.end() && (*p)->GetPriority() <= aProcess->GetPriority();
			p++)
			;
		processTable.insert(p, aProcess);
		MERLIN_PROFILE_ADD_PROCESS(aProcess->GetID());
	}
}

bool ProcessStepManager::RemoveProcess(BunchProcess* aProcess)
{
	proc_itor p = find(processTable.begin(), processTable.end(), aProcess);
	processTable.erase(p);
	return p != processTable.end();
}

void ProcessStepManager::ClearProcesses()
{
	for_each(processTable.begin(), processTable.end(), deleter<BunchProcess>());
	processTable.clear();
}

void ProcessStepManager::SetLogStream(ostream* os)
{
	log = os;
}
