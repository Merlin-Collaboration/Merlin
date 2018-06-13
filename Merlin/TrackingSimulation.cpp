/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <algorithm>
#include <cassert>
#include "MerlinIO.h"
#include "TrackingSimulation.h"

namespace
{

template<class II>
void PerformTracking(ProcessStepManager& aStepper, Bunch& aBunch, bool includeX, bool injOnAxis,
	SimulationOutput* simop, II first, II last)
{
	bool fb = true;
	do
	{
		ComponentFrame* frame = *first;

		if(fb && injOnAxis)
		{
			std::cout << "ignoring first frame transformation" << std::endl;
		}

		if(includeX && !(fb && injOnAxis))
		{
			aBunch.ApplyTransformation(frame->GetEntrancePlaneTransform());
		}

		if(const Transform3D* t = frame->GetEntranceGeometryPatch())
		{
			aBunch.ApplyTransformation(*t);
		}

		if(frame->IsComponent())
		{
			aStepper.Track(frame->GetComponent());
		}
		if(const Transform3D* t = frame->GetExitGeometryPatch())
		{
			aBunch.ApplyTransformation(*t);
		}

		if(includeX)
		{
			aBunch.ApplyTransformation(frame->GetExitPlaneTransform());
		}
		if(simop)
		{
			simop->DoRecord(frame, &aBunch);
		}

		fb = false;
	} while(++first != last);
}

} // end of anonymous namespace

template<class I>
TrackingSimulation::TStepper<I>::TStepper(I f, I l) :
	first(f), last(l), curr(f)
{
	cFrame = *curr;
}

template<class I>
ComponentFrame* TrackingSimulation::TStepper<I>::NextFrame()
{
	curr++;
	return cFrame = curr == last ? nullptr : *curr;
}

TrackingSimulation::TrackingSimulation(const AcceleratorModel::Beamline& bline) :
	bunch(nullptr), incX(true), injOnAxis(false), log(nullptr), handle_me(false), type(beamline), ibunchCtor(nullptr),
	stepper(), theRing(), theBeamline(bline), cstepper(nullptr), simOp(nullptr)
{
}

TrackingSimulation::TrackingSimulation(const AcceleratorModel::RingIterator& aRing) :
	bunch(nullptr), incX(true), injOnAxis(false), log(nullptr), handle_me(false), type(ring), ibunchCtor(nullptr),
	stepper(), theRing(aRing), theBeamline(), cstepper(nullptr), simOp(nullptr)
{
}

TrackingSimulation::TrackingSimulation() :
	bunch(nullptr), incX(true), injOnAxis(false), log(nullptr), handle_me(false), type(undefined), ibunchCtor(nullptr),
	stepper(), theRing(), theBeamline(), cstepper(nullptr), simOp(nullptr)
{
}

void TrackingSimulation::SetBeamline(const AcceleratorModel::Beamline& bline)
{
	theBeamline = bline;
	type = beamline;
}

void TrackingSimulation::SetRing(const AcceleratorModel::RingIterator& aRing)
{
	theRing = aRing;
	type = ring;
}

TrackingSimulation::~TrackingSimulation()
{
	if(bunch)
	{
		delete bunch;
	}
	if(ibunchCtor)
	{
		delete ibunchCtor;
	}
}

Bunch& TrackingSimulation::DoRun(bool new_bunch, bool do_init)
{
	// force initialisation if there is no bunch
	do_init = new_bunch || do_init || (bunch == nullptr);

	if(new_bunch)
	{
		assert(ibunchCtor != nullptr);
		if(bunch != nullptr)
		{
			delete bunch;
		}
		bunch = ibunchCtor->ConstructBunch();
	}

	if(simOp)
	{
		simOp->DoRecordInitialBunch(bunch);
	}

	// NOTE potential bug: if injOnAxis is true then do_init should also be true.
	// do_init == false indicates a continuation of a previous tracking run (a ring)
	// in which case injOnAxis==true may be an error.

	if(injOnAxis && !do_init)
	{
		std::cerr << "*** WARNING: possible tracking error - injOnAxis==true for continued tracking" << std::endl;
	}

	try
	{
		if(do_init)
		{
			stepper.Initialise(*bunch);
		}

		if(type == beamline)
		{
			PerformTracking(stepper, *bunch, incX, injOnAxis, simOp, theBeamline.begin(), theBeamline.end());
		}
		else
		{
			PerformTracking(stepper, *bunch, incX, injOnAxis, simOp, theRing, theRing);
		}
	}
	catch(MerlinException& me)
	{
		if(handle_me)
		{
			MERLIN_ERR << std::endl << me.Msg() << std::endl;
		}
		else
		{
			throw;
		}
	}

	if(simOp)
	{
		simOp->DoRecordFinalBunch(bunch);
	}

	return *bunch;
}

Bunch& TrackingSimulation::Run()
{
	return DoRun(true, true);
}

Bunch& TrackingSimulation::Continue()
{
	return DoRun(false, false);
}

void TrackingSimulation::AddProcess(BunchProcess* proc)
{
	stepper.AddProcess(proc);
}

bool TrackingSimulation::RemoveProcess(BunchProcess* proc)
{
	return stepper.RemoveProcess(proc);
}

void TrackingSimulation::AssumeFlatLattice(bool flat)
{
	incX = !flat;
}

void TrackingSimulation::SetInitialBunchCtor(BunchConstructor* bctor)
{
	//	assert(bctor!=0);
	if(ibunchCtor != nullptr)
	{
		delete ibunchCtor;
	}
	ibunchCtor = bctor;
}

void TrackingSimulation::InitStepper(bool genNewBunch)
{
	assert(ibunchCtor != nullptr);

	if(genNewBunch || bunch == nullptr)
	{
		if(bunch != nullptr)
		{
			delete bunch;
		}
		bunch = ibunchCtor->ConstructBunch();
	}

	stepper.Initialise(*bunch);

	if(cstepper != nullptr)
	{
		delete cstepper;
	}

	if(type == beamline)
	{
		cstepper = new BeamlineStepper(theBeamline.begin(), theBeamline.end());
	}
	else
	//		cstepper = new RingStepper(theRing,theRing);
	{
		cstepper = nullptr;
	}
}

void TrackingSimulation::InitStepper(Bunch* aBunch)
{
	if(bunch != nullptr)
	{
		delete bunch;
	}

	bunch = aBunch;
	stepper.Initialise(*bunch);

	if(cstepper != nullptr)
	{
		delete cstepper;
	}

	if(type == beamline)
	{
		cstepper = new BeamlineStepper(theBeamline.begin(), theBeamline.end());
	}
	else
	{
		cstepper = nullptr;
	}
}

bool TrackingSimulation::StepComponent()
{
	assert(cstepper && cstepper->cFrame);
	ComponentFrame* frame = cstepper->cFrame;
	if(incX)
	{
		bunch->ApplyTransformation(frame->GetEntrancePlaneTransform());
	}
	stepper.Track(frame->GetComponent());
	if(incX)
	{
		bunch->ApplyTransformation(frame->GetExitPlaneTransform());
	}

	if(simOp)
	{
		simOp->DoRecord(frame, bunch);
	}

	return cstepper->NextFrame() != nullptr;
}

void TrackingSimulation::SetOutput(SimulationOutput *simout)
{
	simOp = simout;
}

// SimulationOutput definitions

bool SimulationOutput::IsMember(const string& key)
{
	for(std::vector<StringPattern>::const_iterator p = ids.begin(); p != ids.end(); p++)
	{
		if(p->Match(key))
		{
			return true;
		}
	}
	return false;
}

void SimulationOutput::AddIdentifier(const string& pattern, size_t nocc)
{
	ids.push_back(pattern);
}

void SimulationOutput::DoRecord(const ComponentFrame* frame, const Bunch* bunch)
{
	if(frame->IsComponent())
	{
		string id = (*frame).GetComponent().GetQualifiedName();

		if(output_all || IsMember(id))
		{
			Record(frame, bunch);
		}
	}
}

void SimulationOutput::DoRecordInitialBunch(const Bunch* bunch)
{
	if(output_initial)
	{
		RecordInitialBunch(bunch);
	}
}

void SimulationOutput::DoRecordFinalBunch(const Bunch* bunch)
{
	if(output_final)
	{
		RecordFinalBunch(bunch);
	}
}

void DoRecordFinalBunch(const Bunch* bunch);
