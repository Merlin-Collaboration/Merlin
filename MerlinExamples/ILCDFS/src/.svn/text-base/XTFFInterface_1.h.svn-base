/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef XTFFInterface_h
#define XTFFInterface_h 1

#include "merlin_config.h"
#include <fstream>
#include <string>
#include <set>
#include "AcceleratorModel/AcceleratorModel.h"
#include "BeamModel/BeamData.h"

class AcceleratorModelConstructor;

//  class XTFFInterface
//  Class used to construct a MERLIN model from a MAD TWISS TAPE
//  output listing. XTFF (eXtended Tape File Format) contains
//  SLAC extentions for cavities and acceleration.

class XTFFInterface_1 {
public:

    // Constructer taking the name of the .xtff file, the total bunch
    // charge (particle per bunch) and an option log file. When a non-zero
    // bunch charge is specified, the constructor uses the ELOSS information
    // for the cavities to calculate the reference (matched) energy for the
    // magnet strengths.
    XTFFInterface_1(const std::string& fname, double nb=0, ostream* log =0);
    ~XTFFInterface_1();

    pair<AcceleratorModel*,BeamData*> Parse();
    pair<AcceleratorModel*,BeamData*> Parse(double P_ref);

    // Construct apertures if flag is true (default)
    void IncludeApertures(bool flag) { incApertures = flag; }

    // Treat MAD type as DRIFT
    void TreatTypeAsDrift(const string&);

    // Construct girders
    void ConstructGirders(bool flg) { girders = flg; }

    // data structure for XTFF data
    struct XTFF_Data;

private:

    void ConstructComponent(XTFF_Data&);
    int ParseHeader();
    pair<AcceleratorModel*,BeamData*> Parse1(int);

    std::set<string> driftTypes;

    std::ifstream ifs;
    AcceleratorModelConstructor* mc;
    BeamData* beam0;
    double nb;
    ostream* logos;
    bool incApertures;

    // used for girder construction
    bool girders;
    bool in_g;
};

#endif

