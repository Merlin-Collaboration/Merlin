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

#ifndef BunchProcess_h
#define BunchProcess_h 1

#include "merlin_config.h"
#include <string>

class AcceleratorComponent;
class Bunch;

using namespace std;

//	A single process which is applied to a Bunch object
//	during tracking. A Process can represent any physical
//	(eg. particle transport) or abstract mechanism (eg.
//	output) which may be applied to a bunch at a specific
//	AcceleratorComponent. Concrete processes may be generic
//	to all Bunch model types , or may be specific to a
//	concrete bunch model (see for example ParticleBunch and
//	ParticleBunchTransortProcess).

class BunchProcess
{
public:

    //	Constructing taking the process ID and its priority
    //	(default =0, high).
    explicit BunchProcess (const string& anID, int aPriority = 0);

    virtual ~BunchProcess ();

    //	Set/Get the priority. Highest priority is 0, followed by
    //	1,2,...etc.
    void SetPriority (int p);
    int GetPriority () const;

    //	Initialise this process with the specified Bunch.
    virtual void InitialiseProcess (Bunch& bunch) = 0;

    //	Sets the current accelerator component. This function
    //	should be called just before tracking of the component
    //	begins. Concrete processes should override this function
    //	to perform component and preocess dependent
    //	initialisation.
    virtual void SetCurrentComponent (AcceleratorComponent& component) {
        currentComponent = &component;
    }

    //	Preform the process for the specified step ds.
    virtual void DoProcess (double ds) = 0;

    //	Returns the current maximum step length for this process.
    virtual double GetMaxAllowedStepSize () const = 0;

    //	Returns true if this processes is active.
    bool IsActive () const;

    const string& GetID () const;

protected:

    bool active;
    AcceleratorComponent* currentComponent;

private:

    int priority;
    string ID;
};

// Template class for defining bunch specific BunchProcess classes

template<class B> class TBunchProc : public BunchProcess {
public:

    explicit TBunchProc(const std::string& anID, int aPriority =0)
            : BunchProcess(anID,aPriority)
    {}

    // Sets the current bunch, if bunch is of type bunch
    virtual void InitialiseProcess (Bunch& bunch) {
        currentBunch = dynamic_cast<B*>(&bunch);
    }

protected:

    B* currentBunch;
};

inline BunchProcess::BunchProcess (const string& anID, int aPriority)
        : active(false),priority(aPriority),ID(anID)
{}

inline BunchProcess::~BunchProcess ()
{}

inline void BunchProcess::SetPriority (int p)
{
    priority = p;
}

inline int BunchProcess::GetPriority () const
{
    return priority;
}

inline bool BunchProcess::IsActive () const
{
    return active;
}

inline const string& BunchProcess::GetID () const
{
    return ID;
}

#endif
