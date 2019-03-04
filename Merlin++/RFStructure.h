/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef RFStructure_h
#define RFStructure_h 1

#include "merlin_config.h"
#include "PhysicalConstants.h"
#include "TemplateComponents.h"
#include "RectangularGeometry.h"
#include "RFAcceleratingField.h"

/**
 *	Takes a field which is derived from RFAcceleratingField.
 */

class RFStructure: public TAccCompGF_NC<RectangularGeometry, RFAcceleratingField>
{
public:

	/**
	 *	Returns the frequency.
	 *	@return Frequency
	 */
	double GetFrequency() const;

	/**
	 *	Access the peak amplitude (gradient).
	 */
	void SetAmplitude(double Epk);
	double GetAmplitude() const;

	/**
	 * Accessing the voltage and phase
	 */
	void SetVoltage(double v);
	void SetPhase(double phase);
	void SetVoltagePhasor(const Complex& z);

	double GetVoltage() const;
	double GetPhase() const;
	Complex GetVoltagePhasor() const;

	/**
	 * Calculate the nominal energy gain
	 */
	double GetBeamVoltage() const;

	/**
	 *	Returns the wavelength of the RF (in meter).
	 *	@return Wavelength of RF (m)
	 */
	double GetWavelength() const;

	/**
	 *	Returns the k value (=2pi/wavelength) for the field.
	 *	@return Wavenumber k (\f$ 2\pi/\lambda \f$) for the field
	 */
	double GetK() const;

	/**
	 *	Modify the frequency
	 */
	void SetFrequency(double f);
	void SetWavelength(double lambda);
	void SetK(double k);

protected:

	/**
	 * Protected constructor prevents instantiation of this class.
	 */
	RFStructure(const string& id, double len, RFAcceleratingField* aField);
};

#endif
