// ConstantStrayFieldProcess.h
//
// Merlin 3.0 example code
// v1.0 Nick Walker 21.10.2004
//-------------------------------------------------------------------------


// An example of a BunchProcess which superimposes the effect of a constant
// 'stray' vertical magnetic field on the tracked bunch.
//
// The effect of the field is modelled as a kick which is applied at the
// exit of every element or every distance maxStep (depending which is shorter).
//
// Note this is rather artificial example, but it does illustrate some 
// important concepts.
//-------------------------------------------------------------------------

#ifndef _h_ConstantStrayFieldProcess
#define _h_ConstantStrayFieldProcess 1

#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"

namespace ParticleTracking {
	
	class ConstantStrayFieldProcess : public ParticleBunchProcess {
	public:
		
		// construction
		ConstantStrayFieldProcess(double maxStep, double By);
		
		// Virtual function override 
		virtual void SetCurrentComponent (AcceleratorComponent& component);

		// The following pure virtual function overrides must be provided.
		virtual void DoProcess (double ds);
		virtual double GetMaxAllowedStepSize () const;
		
	private:
		
		double By;
		double maxStep;

		// implementation 
		double kick_ds;
		double s_int;
		double cL;

		void ApplyKick();
	};
	
}; // end namespace ParticleTracking

#endif // _h_ConstantStrayFieldProcess