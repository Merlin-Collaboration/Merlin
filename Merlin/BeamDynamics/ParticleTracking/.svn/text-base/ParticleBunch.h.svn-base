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
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef ParticleBunch_h
#define ParticleBunch_h 1

#include "merlin_config.h"
// PSTypes
#include "BeamModel/PSTypes.h"
// Bunch
#include "BeamModel/Bunch.h"

class Transform3D;

namespace ParticleTracking {

//	Representation of a particle.
typedef PSvector Particle;

//	A Bunch which is represented by an ensemble of
//	(macro-)particles.
class ParticleBunch : public Bunch
{
public:

    typedef Particle particle_type;
    typedef PSvectorArray::iterator iterator;
    typedef PSvectorArray::const_iterator const_iterator;

    //	Constructs a ParticleBunch using the specified momentum,
    //	total charge and the particle array. Note that on exit,
    //	particles is empty.
    ParticleBunch (double P0, double Q, PSvectorArray& particles);

    //	Read phase space vectors from specified input stream.
    ParticleBunch (double P0, double Q, std::istream& is);

    //	Constructs an empty ParticleBunch with the specified
    //	momentum P0 and charge per macro particle Qm (default =
    //	+1).
    ParticleBunch (double P0, double Qm = 1);

    //	Returns the total charge (in units of e).
    virtual double GetTotalCharge () const;

    virtual PSmoments& GetMoments (PSmoments& sigma) const;
    virtual PSmoments2D& GetProjectedMoments (PScoord u, PScoord v, PSmoments2D& sigma) const;
    virtual PSvector& GetCentroid (PSvector& p) const;
    virtual Point2D GetProjectedCentroid (PScoord u, PScoord v) const;

    // Calculate mean (first) and rms (second) of the specified coordinate.
    std::pair<double,double> GetMoments(PScoord u) const;

    //	Set the reference momentum to the mean (centroid)
    //	momentum of the bunch. Returns the new value in GeV/c.
    virtual double AdjustRefMomentumToMean ();

    //	Adjusts the reference moment by the relative dpp.
    virtual double AdjustRefMomentum (double dpp);

    //	Set the reference time to the mean (centroid) time of
    //	the bunch. Returns the new value in ct (meters).
    virtual double AdjustRefTimeToMean ();

    //	Used to generate a 1-D profile of the bunch projected
    //	onto the specified coordinate. The total area of the
    //	historgram is normalised to unity.
    virtual Histogram& ProjectDistribution (PScoord axis, Histogram& hist) const;

    //	Apply the specified 3D coordinate transformation to the
    //	bunch. Returns true if successful (note that it may not
    //	be possible to apply a general 3D transformation to some
    //	concrete bunch representations).
    virtual bool ApplyTransformation (const Transform3D& t);

    //	Sorts the particles into ascending ct order.
    virtual void SortByCT ();

    //	Output a bunch-model dependent representation to the
    //	specified stream.
    virtual void Output (std::ostream& os) const;

    //	Add a (macro-)particle to the bunch.
    virtual size_t AddParticle (const Particle& p);

    //	Sets the particle charge per macro-particle.
    void SetMacroParticleCharge (double q);

    ParticleBunch::iterator begin ();
    ParticleBunch::iterator end ();
    ParticleBunch::const_iterator begin () const;
    ParticleBunch::const_iterator end () const;
    size_t size () const;
    virtual void push_back (const Particle& p);
    virtual ParticleBunch::iterator erase (ParticleBunch::iterator p);

    PSvectorArray& GetParticles ();
    const PSvectorArray& GetParticles () const;

    //	Returns the first particle in the bunch.
    const Particle& FirstParticle () const;
    Particle& FirstParticle ();

    //	Sets the centroid of the particle bunch to be exactly x0.
    void SetCentroid (const Particle& x0);

protected:

    PSvectorArray pArray;

private:

    //	Charge per macro-particle
    double qPerMP;
};

inline size_t ParticleBunch::AddParticle (const Particle& p)
{
    pArray.push_back(p);
    return size();
}

inline void ParticleBunch::SetMacroParticleCharge (double q)
{
    qPerMP = q;
    SetChargeSign(q);
}

inline ParticleBunch::iterator ParticleBunch::begin ()
{
    return pArray.begin();
}

inline ParticleBunch::iterator ParticleBunch::end ()
{
    return pArray.end();
}

inline void ParticleBunch::push_back (const Particle& p)
{
    AddParticle(p);
}

inline ParticleBunch::const_iterator ParticleBunch::begin () const
{
    return pArray.begin();
}

inline ParticleBunch::const_iterator ParticleBunch::end () const
{
    return pArray.end();
}

inline size_t ParticleBunch::size () const
{
    return pArray.size();
}

inline ParticleBunch::iterator ParticleBunch::erase (ParticleBunch::iterator p)
{
    return pArray.erase(p);
}

inline PSvectorArray& ParticleBunch::GetParticles ()
{
    return pArray;
}

inline const PSvectorArray& ParticleBunch::GetParticles () const
{
    return pArray;
}

inline const Particle& ParticleBunch::FirstParticle () const
{
    return pArray.front();
}

inline Particle& ParticleBunch::FirstParticle ()
{
    return pArray.front();
}

}; // end namespace ParticleTracking

#endif
