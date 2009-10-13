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

#ifndef ProcessStepManager_h
#define ProcessStepManager_h 1

#include "merlin_config.h"
#include <list>
#include <ostream>

class AcceleratorComponent;
class BunchProcess;
class Bunch;

//	Responsible for coordinating the tracking of a bunch
//	through a single AcceleratorComponent. Tracking occurs
//	by the application of a series of processes, which
//	effectively intergrate the bunch motion along the
//	beamline.

class ProcessStepManager
{
public:

    //	Construction/destruction.
    ProcessStepManager ();
    ~ProcessStepManager ();

    //	Initialise the step manager with the specified (initial)
    //	bunch.
    void Initialise (Bunch& bunch);

    //	Track the specified component. The current bunch object
    //	is updated accordingly.
    void Track (AcceleratorComponent& component);

    //	Returns the total length integrated since the last call
    //	to Initialise(Bunch&).
    double GetIntegratedLength ();

    //	Add a process.
    void AddProcess (BunchProcess* aProcess);

    //	Remove aProcess from the current process table. Returns
    //	true if aProcess was present, otherwise false.
    bool RemoveProcess (BunchProcess* aProcess);

    //	Remove and destroys all processes in the current process
    //	table.
    void ClearProcesses ();

    void SetLogStream (std::ostream* os);

private:

    //	The current intgrated length.
    double total_s;
    std::ostream* log;
    //	list of processes in order of priority.
    //	ordered
    std::list<BunchProcess*> processTable;
};

#endif
