/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2007/03/23 14:58:51 $
// $Revision: 1.8 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef TrackingSimulation_h
#define TrackingSimulation_h 1

#include "merlin_config.h"
// Bunch
#include "BeamModel/Bunch.h"
// BunchConstructor
#include "BeamModel/BunchConstructor.h"
// ProcessStepManager
#include "BeamDynamics/ProcessStepManager.h"
// BunchProcess
#include "BeamDynamics/BunchProcess.h"
// AcceleratorModel
#include "AcceleratorModel/AcceleratorModel.h"
// ComponentFrame
#include "AcceleratorModel/Frames/ComponentFrame.h"
// MerlinException
#include "Exception/MerlinException.h"

// Abstract class used to record bunch and frame information
// after every component frame tracked

class SimulationOutput {
public:
    virtual ~SimulationOutput() {}

    void DoRecord(const ComponentFrame* frame, const Bunch* bunch);
    void DoRecordInitialBunch(const Bunch* bunch);
    void DoRecordFinalBunch(const Bunch* bunch);

	// Output control
	void AddIdentifier(const std::string& pattern, size_t nocc =1);

	// Public output flags
	bool output_all;
	bool output_initial;
	bool output_final;

protected:

	virtual void Record(const ComponentFrame* frame, const Bunch* bunch) =0;
    virtual void RecordInitialBunch(const Bunch* bunch)=0;
    virtual void RecordFinalBunch(const Bunch* bunch)=0;

private:

	bool IsMember(const std::string&);
	std::vector<StringPattern> ids;
};

//	A beam dynamics simulation. TrackingSimulation tracks a
//	Bunch object through a specified Beamline (single pass).
//	Tracking is performed by one or more BunchProcess
//	objects, one of which is generally a transport process
//	which is responsible for "tracking" the concrete Bunch
//	implementation through the Beamline components. The
//	initial Bunch is generated on each Run using the
//	supplied BunchConstructor.

class TrackingSimulation
{
private:

    typedef enum {beamline, ring, undefined} lattice_type;

    //	Base class for stepping (iterating) through a beamline
    //	or ring.
    class Stepper
    {
    public:
        virtual ComponentFrame* NextFrame () = 0;
        ComponentFrame* cFrame;
    };

    //	Stepper implementation to take either a ring or a
    //	beamline iterator.
    template <class I>
    class TStepper : public Stepper
    {
    public:
        TStepper (I f, I l);
        virtual ComponentFrame* NextFrame ();
    private:
        I first;
        I last;
        I curr;
    };

public:

    typedef TStepper< AcceleratorModel::BeamlineIterator  > BeamlineStepper;
    typedef TStepper< AcceleratorModel::RingIterator  > RingStepper;

    // Construct a TrackingSimulation without specifying a beamline
    TrackingSimulation();
    TrackingSimulation (const AcceleratorModel::Beamline& bline);
    TrackingSimulation (const AcceleratorModel::RingIterator& aRing);

    virtual ~TrackingSimulation ();

    // Set the beamline or ring to be used by the tracker
    void SetBeamline(const AcceleratorModel::Beamline& bline);
    void SetRing(const AcceleratorModel::RingIterator& aRing);

    //	Run the simulation. A new bunch is constructed using
    //  the current BunchConstructor object, and all the 
    //  processes are initialised.
    virtual Bunch& Run();

    //  Run the simulation again without constructing a new bunch
    //  or (re-)initialising the BunchProcesses. This function can
    //  only be called after an initial Run() method is called.
    virtual Bunch& Continue();

    //	Add a BunchProcess to the simulation.
    void AddProcess (BunchProcess* proc);

    //	Remove the specified BunchProcess from the simulation.
    //	Returns true if successful.
    bool RemoveProcess (BunchProcess* proc);

    //	if flat is true, the simulation ignores all coordinate
    //	transformations between lattice components. If flat is
    //	false (default), the simulation applies all necessary
    //	in- and out-transformations for each component in the
    //	beamline.
    void AssumeFlatLattice (bool flat);

	//  If onAxis is true, tracking simulation ignores any coordination 
	//  transformation for the first component frame tracked. This effectively
	//  injects the 'beam' on the local components axis, irrespective of its 
	//  alignment state.
	void InjectBeamOnAxis(bool onAxis) {
		injOnAxis = onAxis;
	}

    //	Sets the initial bunch constructor.
    void SetInitialBunchCtor (BunchConstructor* bctor);

    void SetLogStream (ostream& os)
    {
        log = &os;
    }

    void SetLoggingOn ()
    {
        stepper.SetLogStream(log);
    }

    void SetLoggingOff ()
    {
        stepper.SetLogStream(0);
    }

    //	If handle==true, then any MerlinException thrown while
    //	the simulation is running is handled locally, and Run()
    //	command exits normally. If false, any exception is
    //	re-thrown.
    void HandleExceptions (bool handle)
    {
        handle_me = handle;
    }

    //	The following methods are for use when stepping though
    //	the beamline.
    //
    //	Initialise a new bunch object for tracking, and resets
    //	the current frame to the first frame of the beamline.
    void InitStepper (bool genNewBunch = true);

    //	Sets the specified bunch for tracking, and resets
    //	the current frame to the first frame of the beamline.
    void InitStepper (Bunch* aBunch);

    //	Track the bunch though the next ComponentFrame in the
    //	beamline. Returns false if that component is the last in
    //	the beamline.
    bool StepComponent ();

    //	Returns the current frame, i.e. that frame that is about
    //	to be tracked by a call to StepComponent().
    ComponentFrame& GetCurrentFrame ()
    {
        assert(cstepper);
        return *(cstepper->cFrame);
    }

    //	Returns the current component, i.e. that component that
    //	is about to be tracked by a call to StepComponent().
    AcceleratorComponent& GetCurrentComponent ()
    {
        assert(cstepper);
        return cstepper->cFrame->GetComponent();
    }

    void SetOutput(SimulationOutput* simout);

protected:

    Bunch* bunch;
    Bunch& DoRun(bool new_bunch, bool do_init);

private:

    bool incX;
	bool injOnAxis;
    ostream* log;
    bool handle_me;
    lattice_type type;

    BunchConstructor* ibunchCtor;
    ProcessStepManager stepper;

    AcceleratorModel::RingIterator theRing;
    AcceleratorModel::Beamline theBeamline;
    Stepper* cstepper;
    SimulationOutput* simOp;
};

//	A template class which is used to implement a bunch
//	implementation specific TrackingSimulation. P must be a
//	bunch specific transport process which has only a
//	default constructor. It must also define P::bunch_type
//	which specifies the type of bunch implementation to be
//	tracked.

template <class P>
class BunchTracker : public TrackingSimulation
{
public:

    typedef typename P::bunch_type bunch_type;
    typedef P transport_process;

    //	Constructor taking the beamline to be tracked. The
    //	initial bunch must be specified by supplying a concrete
    //	BunchConstructor via SetInitialBunchCtor.
    explicit BunchTracker (const AcceleratorModel::Beamline& bline)
            : TrackingSimulation(bline)
    {
        AddProcess(new P());
    }

    //	Constructs a BunchTracker taking the beamline to be
    //	tracked and the initial bunch to be used.
    BunchTracker (const AcceleratorModel::Beamline& bline, bunch_type* bunch0)
            : TrackingSimulation(bline)
    {
        AddProcess(new P());
        SetInitialBunchCtor(new StaticBunchCtor<bunch_type>(bunch0));
    }

    bunch_type& GetBunch ()
    {
        return static_cast<bunch_type&>(*bunch);
    }

    transport_process& GetTransportProcess ()
    {
        return *transport;
    }

private:

    transport_process* transport;
};


#endif
