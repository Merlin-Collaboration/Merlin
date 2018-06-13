/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_TTrackSim
#define _h_TTrackSim 1

#include "TrackingSimulation.h"
#include "ComponentTracker.h"

/**
 * Transport process template function.
 * Defines a transport BunchProcess which uses a concrete ComponentTracker
 * class to track provide the tracking.
 *
 * Template argument T must be a valid concrete ComponentTracker class.
 */
template<class T>
class TTrnsProc: public TBunchProc<__TYPENAME__ T::bunch_type>
{
public:

	typedef __TYPENAME__ T::bunch_type bunch_type;
	typedef __TYPENAME__ T::B_Integrator integrator_type;
	typedef __TYPENAME__ T::ISetBase integrator_set_base;

	/**
	 * Construct process using the default integrator set
	 */
	TTrnsProc() :
		TBunchProc<__TYPENAME__ T::bunch_type>("TRANSPORT", 0), ctracker()
	{
	}

	/**
	 * Construct process using an explicit integrator set
	 */
	TTrnsProc(const integrator_set_base& iset) :
		TBunchProc<__TYPENAME__ T::bunch_type>("TRANSPORT", 0), ctracker(iset)
	{
	}

	/**
	 * Register an addition (or override) integrator.
	 */
	bool RegisterIntegrator(integrator_type* intg)
	{
		return ctracker.Register(intg);
	}

	void InitialiseProcess(Bunch& bunch)
	{
		TBunchProc<bunch_type>::InitialiseProcess(bunch);
		if(this->currentBunch != nullptr)
		{
			this->active = true;
			ctracker.SetBunch(*this->currentBunch);
		}
	}

	void SetCurrentComponent(AcceleratorComponent& component)
	{
		if(this->active)
		{
			component.PrepareTracker(ctracker);
			this->currentComponent = &component;
		}
		else
		{
			this->currentComponent = nullptr;
		}
	}

	void DoProcess(double ds)
	{
		if(this->active)
		{
			ctracker.TrackStep(ds);
		}
	}

	double GetMaxAllowedStepSize() const
	{
		return this->active ? ctracker.GetRemainingLength() : 0;
	}

	void SetIntegratorSet(const integrator_set_base* iset)
	{
		ctracker.ClearIntegratorSet();
		iset->Init(ctracker);
	}

private:

	T ctracker;
};

/**
 * template class to instantiate a TrackingSimulation class which
 * uses a specific concrete ComponentTracker to track the
 * associated bunch representation through a beamline.
 *
 * Template parameter T must be a valid concrete ComponentTracker class
 */
template<class T>
class TTrackSim: public TrackingSimulation
{
public:

	typedef __TYPENAME__ T::bunch_type bunch_type;
	typedef __TYPENAME__ bunch_type::particle_type particle_type;
	typedef TTrnsProc<T> transport_process;
	typedef __TYPENAME__ transport_process::integrator_type integrator_type;
	typedef __TYPENAME__ T::ISetBase integrator_set_base;

	/**
	 * Constructor taking the beamline to be tracked and a
	 * pointer to the initial ParticleBunch. If bunch0=0
	 * (default) the initial beam must be specified by either a
	 * call to SetInitialBeam(), or to SetInitialBeamCtor (from
	 * TrackingSimulation).
	 */
	explicit TTrackSim(const AcceleratorModel::Beamline& bline, bunch_type* bunch0 = nullptr, bool del = false);

	/**
	 * Constructor used for single particle tracking. The
	 * constructor takes the beamline to be tracked, and the
	 * initial particle.
	 */
	TTrackSim(const AcceleratorModel::Beamline& bline, const particle_type& p, double Pref);

	/**
	 * Constructor taking the beamline to be tracked and a
	 * pointer to the initial ParticleBunch. If bunch0=0
	 * (default) the initial beam must be specified by either a
	 * call to SetInitialBeam(), or to SetInitialBeamCtor (from
	 * TrackingSimulation).
	 */
	explicit TTrackSim(const AcceleratorModel::RingIterator& ring, bunch_type* bunch0 = nullptr, bool del = false);

	/**
	 * Constructor used for single particle tracking. The
	 * constructor takes the beamline to be tracked, and the
	 * initial particle.
	 */
	TTrackSim(const AcceleratorModel::RingIterator& ring, const particle_type& p, double Pref);

	/**
	 * Default constructor
	 */
	TTrackSim();

	/**
	 * Register an additional (or override) integrator.
	 */
	bool RegisterIntegrator(integrator_type* intg)
	{
		return transportProc->RegisterIntegrator(intg);
	}

	/**
	 * Sets the initial ParticleBunch for future tracking operations.
	 */
	void SetInitialBunch(bunch_type* pbunch0, bool del = false)
	{
		if(pbunch0 != nullptr)
		{
			SetInitialBunchCtor(MakeBunchCtor(pbunch0, del));
		}
	}

	/**
	 * Overrides the current bunch constructor and
	 * tracks the supplied bunch.
	 */
	bunch_type* Track(bunch_type*);

	/**
	 * Sets the initial particle for single-particle tracking.
	 */
	void SetInitialParticle(const particle_type& p, double Pref);

	/**
	 * Returns a reference to the current tracked bunch.
	 */
	const bunch_type& GetTrackedBunch() const
	{
		return static_cast<const bunch_type&>(*bunch);
	}

	bunch_type& GetTrackedBunch()
	{
		return static_cast<bunch_type&>(*bunch);
	}

	void SetIntegratorSet(const integrator_set_base* iset)
	{
		transportProc->SetIntegratorSet(iset);
	}

private:
	transport_process* transportProc;

	//Copy protection
	TTrackSim(const TTrackSim& rhs);
	TTrackSim& operator=(const TTrackSim& rhs);

};

// template implementation
/**
 * Standard tracker constructor
 */
template<class T>
TTrackSim<T>::TTrackSim(const AcceleratorModel::Beamline& bline, bunch_type* bunch0, bool del) :
	TrackingSimulation(bline), transportProc(new transport_process())
{
	SetInitialBunch(bunch0, del);
	AddProcess(transportProc);
}

template<class T>
TTrackSim<T>::TTrackSim(const AcceleratorModel::Beamline& bline, const particle_type& p, double Pref) :
	TrackingSimulation(bline), transportProc(new transport_process())
{
	SetInitialParticle(p, Pref);
	AddProcess(transportProc);
}

template<class T>
TTrackSim<T>::TTrackSim(const AcceleratorModel::RingIterator& ring, bunch_type* bunch0, bool del) :
	TrackingSimulation(ring), transportProc(new transport_process())
{
	SetInitialBunch(bunch0, del);
	AddProcess(transportProc);
}

template<class T>
TTrackSim<T>::TTrackSim(const AcceleratorModel::RingIterator& ring, const particle_type& p, double Pref) :
	TrackingSimulation(ring), transportProc(new transport_process())
{
	SetInitialParticle(p, Pref);
	AddProcess(transportProc);
}

template<class T>
TTrackSim<T>::TTrackSim() :
	TrackingSimulation(), transportProc(new transport_process())
{
	AddProcess(transportProc);
}

template<class T>
void TTrackSim<T>::SetInitialParticle(const particle_type& p, double Pref)
{
	bunch_type* b = new bunch_type(Pref);
	b->AddParticle(p);
	SetInitialBunch(b, true);
}
/**
 * Standard Track "Tracker"
 */
template<class T>
__TYPENAME__ TTrackSim<T>::bunch_type * TTrackSim<T>::Track(bunch_type * aBunch)
{
	if(bunch != nullptr)
	{
		delete bunch;
	}

	bunch = aBunch;
	try
	{
		DoRun(false, true);
	}
	catch(...)
	{
		bunch = nullptr;
		throw;
	}

	bunch = nullptr;
	return aBunch;
}

#endif
