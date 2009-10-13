//## begin module%1.4%.codegen_version preserve=yes
//   Read the documentation to learn more about C++ code generator
//   versioning.
//## end module%1.4%.codegen_version

//## begin module%3D206832012A.cm preserve=no
/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file Merlin\BeamDynamics\Utilities\ComponentStepper.cpp
 * last modified 01/07/02 16:36:22
 */
//## end module%3D206832012A.cm

//## begin module%3D206832012A.cp preserve=no
/*
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for 
 * Charge Particle Accelerator Simulations
 * Copyright (c) 2001 by The Merlin Collaboration.
 * - ALL RIGHTS RESERVED - 
 *
 * Permission to use, copy, modify, distribute and sell this
 * software and its documentation for any purpose is hereby
 * granted without fee, provided that the above copyright notice
 * appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation.
 * No representations about the suitability of this software for
 * any purpose is made. It is provided "as is" without express
 * or implied warranty.
 */
//## end module%3D206832012A.cp

//## Module: ComponentStepper%3D206832012A; Package body
//## Subsystem: Merlin::BeamDynamics::Utilities%3D2067F3016F
//## Source file: C:\C++\Merlin\BeamDynamics\Utilities\ComponentStepper.cpp

//## begin module%3D206832012A.includes preserve=yes
#include "NumericalUtils/utils.h"
//## end module%3D206832012A.includes

// AcceleratorComponent
#include "AcceleratorModel/AcceleratorComponent.h"
// ComponentStepper
#include "BeamDynamics/Utilities/ComponentStepper.h"
//## begin module%3D206832012A.additionalDeclarations preserve=yes
//## end module%3D206832012A.additionalDeclarations


// Class ComponentDivider 

//## Operation: ComponentDivider%3D20632A03C7
ComponentDivider::ComponentDivider (int ns, double min_step)
  //## begin ComponentDivider::ComponentDivider%3D20632A03C7.initialization preserve=yes
  :s(0),next_s(0),delta_s(0),minStep(min_step),nstep(ns)
  //## end ComponentDivider::ComponentDivider%3D20632A03C7.initialization
{
  //## begin ComponentDivider::ComponentDivider%3D20632A03C7.body preserve=yes
  //## end ComponentDivider::ComponentDivider%3D20632A03C7.body
}



//## Other Operations (implementation)
//## Operation: SetComponent%3D20634C03E4
void ComponentDivider::SetComponent (AcceleratorComponent& cmp)
{
  //## begin ComponentDivider::SetComponent%3D20634C03E4.body preserve=yes
	double l = cmp.GetLength();
	if(l!=0) {
		delta_s = l/nstep;
		if(delta_s<minStep) {
			int ns = (l/minStep)+1;
			delta_s = l/ns;
		}
		next_s = delta_s;
	}
	else
		next_s = 0;
	s=0;
  //## end ComponentDivider::SetComponent%3D20634C03E4.body
}

//## Operation: Increment%3D20634D0006
bool ComponentDivider::Increment (double ds)
{
  //## begin ComponentDivider::Increment%3D20634D0006.body preserve=yes
	s+=ds;
	if(fequal(s,next_s)) {
		next_s+=delta_s;
		return true;
	}
	return false;
  //## end ComponentDivider::Increment%3D20634D0006.body
}

//## Operation: DistanceToStepBoundary%3D20634D0010
double ComponentDivider::DistanceToStepBoundary () const
{
  //## begin ComponentDivider::DistanceToStepBoundary%3D20634D0010.body preserve=yes
	return next_s-s;
  //## end ComponentDivider::DistanceToStepBoundary%3D20634D0010.body
}

//## begin module%3D206832012A.epilog preserve=yes
//## end module%3D206832012A.epilog
