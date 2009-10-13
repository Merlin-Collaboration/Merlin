/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:53 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef h_SliceMPBunch
#define h_SliceMPBunch 1

#include "merlin_config.h"
#include <iostream>

// ReferenceParticle
#include "BeamModel/Bunch.h"
#include "BeamDynamics/SMPTracking/SliceMacroParticle.h"
#include <vector>

class Transform3D;
class Histogram;

// A bunch which is represented by a collection of
// sliced macro-particles (SliceMacroParticle) as used
// in typical linac beam dynamics simulations.
//

namespace SMPTracking {

class SMPBunch : public  Bunch  {
public:

    typedef SliceMacroParticle particle_type;
    typedef std::vector<SliceMacroParticle> SliceMPArray;
    typedef SliceMPArray::iterator iterator;
    typedef SliceMPArray::const_iterator const_iterator;
    typedef particle_type value_type; // needed for mapping templates

    //	Constructor taking the initial reference momentum in
    //	GeV/c, and the total charge/e of the bunch.
    SMPBunch (double p, double q);

    //    Construct the bunch from the specified file. The
    //    file format is compatible with that used by the Output
    //    method.
    SMPBunch (const std::string& fname);

    //	virtual destructor.
    virtual ~SMPBunch ();

    //	Returns the total charge (in units of e).
    virtual double GetTotalCharge () const;

    //    Calculation of first- and second-order moments
    virtual PSmoments& GetMoments (PSmoments& sigma) const;
    virtual PSmoments2D& GetProjectedMoments (PScoord u, PScoord v, PSmoments2D& sigma) const;
    virtual PSvector& GetCentroid (PSvector& p) const;
    virtual Point2D GetProjectedCentroid (PScoord u, PScoord v) const;

    // increment the centroid of the bunch by the specified amount
    void AdjustCentroid(const PSvector&);

    //	Set the reference momentum to the mean (centroid)
    //	momentum of the bunch. Returns the new value in GeV/c.
    virtual double AdjustRefMomentumToMean ();

    //	Set the reference time to the mean (centroid) time of
    //	the bunch. Returns the new value in ct (meters).
    virtual double AdjustRefTimeToMean ();

    //	Output the slices to the stream. The flat ascii format
    //    is one row per SliceMacroParticle with the following
    //    format:
    //      row 1 cols 1-3: P0[GeV/c] z0[m] Np
    //      rows 2-(Np+1)
    //      columns 1,3  : q[e] z[m] dp/p
    //      columns 4-7  : mean values for x x' y y' in [m] and [r]
    //      columns 8-11 : RMS values for  x x' y y' in [m] and [r]
    //      columns 12-17: correlations <xx'> <yx> <yx'> <y'x> <y'x'> <y'y>
    virtual void Output (std::ostream& os) const;

    //	Used to generate a 1-D profile of the bunch projected
    //	onto the specified coordinate. The total area of the
    //	historgram is normalised to unity.
    virtual Histogram& ProjectDistribution (PScoord axis, Histogram& hist) const;

    //	Apply the specified 3D coordinate transformation to the
    //	bunch. Note that SliceMacroParticles do not support general
    //    transformations, only a sub-set. Returns true if the
    //    the transformation is supported, otherwise false.
    virtual bool ApplyTransformation (const Transform3D& t);

    // Add a macro-particles to the bunch. do_sort==true causes a
    // call to SortByCT().
    void AddParticle(const SliceMacroParticle& p, bool do_sort=false);

    // Sorts the macro-particles in ascending ct order.
    void SortByCT();

    //    Accessors and iterators
    const SliceMPArray& GetSlices() const { return slices; }
    size_t Size() const { return slices.size(); }
    SliceMacroParticle& Get(size_t i) { return slices[i]; }
    const SliceMacroParticle& Get(size_t i) const { return slices[i]; }

    iterator begin() { return slices.begin(); }
    iterator end() { return slices.end(); }
    const_iterator begin() const { return slices.begin(); }
    const_iterator end() const { return slices.end(); }

private:

    double Qt;
    SliceMPArray slices;
};

}; // end namespace SMPTracking

#endif
