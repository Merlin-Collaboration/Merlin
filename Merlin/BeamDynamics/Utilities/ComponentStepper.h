//## begin module%1.4%.codegen_version preserve=yes
//   Read the documentation to learn more about C++ code generator
//   versioning.
//## end module%1.4%.codegen_version

//## begin module%3D2067E9003E.cm preserve=no
/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (1999)
 * 
 * file Merlin\BeamDynamics\Utilities\ComponentStepper.h
 * last modified 01/07/02 16:36:22
 */
//## end module%3D2067E9003E.cm

//## begin module%3D2067E9003E.cp preserve=no
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
//## end module%3D2067E9003E.cp

//## Module: ComponentStepper%3D2067E9003E; Package specification
//## Subsystem: Merlin::BeamDynamics::Utilities%3D2067F3016F
//## Source file: C:\C++\Merlin\BeamDynamics\Utilities\ComponentStepper.h

#ifndef ComponentStepper_h
#define ComponentStepper_h 1

//## begin module%3D2067E9003E.additionalIncludes preserve=no
#include "merlin_config.h"
//## end module%3D2067E9003E.additionalIncludes

//## begin module%3D2067E9003E.includes preserve=yes
//## end module%3D2067E9003E.includes


class AcceleratorComponent;

//## begin module%3D2067E9003E.additionalDeclarations preserve=yes
//## end module%3D2067E9003E.additionalDeclarations


//## Class: ComponentStepper%3D205B1202F5; Abstract
//	Utility class used to calculate the required steps
//	through a specific component.
//## Category: Merlin::BeamDynamics::Utilities%3D205B070385
//## Subsystem: Merlin::BeamDynamics::Utilities%3D2067F3016F
//## Persistence: Transient
//## Cardinality/Multiplicity: n



//## Uses: <unnamed>%3D2068A801B6;AcceleratorComponent { -> F}

class ComponentStepper 
{
  public:

    //## Other Operations (specified)
      //## Operation: SetComponent%3D205B5500FD
      virtual void SetComponent (AcceleratorComponent& cmp) = 0;

      //## Operation: Increment%3D205B730092
      //	Increments the step distance and returns true if on a
      //	step boundary, otherwise false.
      virtual bool Increment (double ds) = 0;

      //## Operation: DistanceToStepBoundary%3D205EE8004E
      //	Returns the distance to the next step boundary.
      virtual double DistanceToStepBoundary () const = 0;

  protected:
  private:
  private: //## implementation
};

//## Class: ComponentDivider%3D205F5500F5
//	A class of stepper which divides a non-zero length
//	component into number of equal steps.
//## Category: Merlin::BeamDynamics::Utilities%3D205B070385
//## Subsystem: Merlin::BeamDynamics::Utilities%3D2067F3016F
//## Persistence: Transient
//## Cardinality/Multiplicity: n



class ComponentDivider : public ComponentStepper  //## Inherits: <unnamed>%3D205FC8037B
{
  public:
    //## Constructors (specified)
      //## Operation: ComponentDivider%3D20632A03C7
      //	Constructor taking the number of steps to take per
      //	component, together with the minimum step distance
      //	(default 0). If for a given component, the calculate
      //	step length is less than min_step, then the step number
      //	is adjusted to give the closest number of integral steps
      //	to match min_step.
      ComponentDivider (int ns, double min_step = 0);


    //## Other Operations (specified)
      //## Operation: SetComponent%3D20634C03E4
      virtual void SetComponent (AcceleratorComponent& cmp);

      //## Operation: Increment%3D20634D0006
      //	Increments the step distance and returns true if on a
      //	step boundary, otherwise false.
      virtual bool Increment (double ds);

      //## Operation: DistanceToStepBoundary%3D20634D0010
      //	Returns the distance to the next step boundary.
      virtual double DistanceToStepBoundary () const;

  protected:
  private:
    // Data Members for Class Attributes

      //## Attribute: s%3D20639200AF
      //	Current component length
      //## begin ComponentDivider::s%3D20639200AF.attr preserve=no  private: double {UA} 
      double s;
      //## end ComponentDivider::s%3D20639200AF.attr

      //## Attribute: next_s%3D206397002A
      //	Integrated length
      //## begin ComponentDivider::next_s%3D206397002A.attr preserve=no  private: double {UA} 
      double next_s;
      //## end ComponentDivider::next_s%3D206397002A.attr

      //## Attribute: delta_s%3D2063CD03D5
      //## begin ComponentDivider::delta_s%3D2063CD03D5.attr preserve=no  private: double {UA} 
      double delta_s;
      //## end ComponentDivider::delta_s%3D2063CD03D5.attr

      //## Attribute: minStep%3D20641902B2
      //## begin ComponentDivider::minStep%3D20641902B2.attr preserve=no  private: double {UA} 
      double minStep;
      //## end ComponentDivider::minStep%3D20641902B2.attr

      //## Attribute: nstep%3D2064240272
      //## begin ComponentDivider::nstep%3D2064240272.attr preserve=no  private: int {UA} 
      int nstep;
      //## end ComponentDivider::nstep%3D2064240272.attr

  private: //## implementation
};

// Class ComponentStepper 

// Class ComponentDivider 

// Class ComponentStepper 

//## begin module%3D2067E9003E.epilog preserve=yes
//## end module%3D2067E9003E.epilog


#endif
