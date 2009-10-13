/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/20 13:42:54 $
// $Revision: 1.6 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef RFStructure_h
#define RFStructure_h 1

#include "merlin_config.h"
#include "NumericalUtils/PhysicalConstants.h"
// TemplateComponents
#include "AcceleratorModel/StdComponent/TemplateComponents.h"
// RectangularGeometry
#include "AcceleratorModel/StdGeometry/RectangularGeometry.h"
// RFAcceleratingField
#include "AcceleratorModel/StdField/RFAcceleratingField.h"

//	Takes a field which is derived from RFAcceleratingField.

class RFStructure : public TAccCompGF_NC<RectangularGeometry,RFAcceleratingField>
{
public:

    //	Returns the frequency.
    double GetFrequency () const;

    //	Access the peak amplitude (gradient).
    void SetAmplitude (double Epk);
    double GetAmplitude () const;

	// Accessing the voltage and phase
	void SetVoltage(double v);
    void SetPhase (double phase);
	void SetVoltagePhasor(const Complex& z);

	double GetVoltage() const;
    double GetPhase () const;
	Complex GetVoltagePhasor() const;

    // Calculate the nominal energy gain
    double GetBeamVoltage() const;

    //	Returns the wavelength of the RF (in meter).
    double GetWavelength () const;

    //	Returns the k value (=2pi/wavelength) for the field.
    double GetK () const;

    //	Modify the frequency
    void SetFrequency (double f);
    void SetWavelength (double lambda);
    void SetK (double k);

protected:

	// Protected constructor prevents insantiation of this class.
    RFStructure (const string& id, double len, RFAcceleratingField* aField);
};

#endif
