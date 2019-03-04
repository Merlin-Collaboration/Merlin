/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Klystron_h
#define Klystron_h 1

#include <string>
#include <vector>
#include <complex>
typedef std::complex<double> Complex;

#include "merlin_config.h"
#include "ModelElement.h"

/**
 * Represents a Klystron, which can be connected to one or
 * more RF cavities. A Klystron has a single Voltage and Phase
 * which is then applied to the attached cavities in one of
 * two modes:
 *
 * vector_sum:
 *    maintains any relative differences in
 *    local phase and voltage of the associated cavities, but
 *    scales them linearly to achieve the total correct vector
 *    sum.
 *
 * balanced:
 *    Sets each associated cavity to the identical phase and
 *    amplitude. If the request phase and amplitude of the klystron
 *    is set to the complex number z, then each of the n cavities
 *    is set to z/n.
 */

class RFStructure;

class Klystron: public ModelElement
{
public:

	enum Mode
	{
		vector_sum,
		balanced

	};

	/**
	 *	Constructor taking the name of the element.
	 */
	Klystron(const std::string& aName, const std::vector<RFStructure*>& cavs, Mode m = balanced);

	virtual ~Klystron();

	// Accessors
	double GetVoltage() const;
	double GetPhase() const;
	Complex GetVoltagePhasor() const;
	void SetVoltage(double v);
	void SetPhase(double phi);
	void SetVoltagePhasor(const Complex&);

	size_t GetNumberOfCavities() const
	{
		return rf_cavs.size();
	}

	/**
	 *	Returns the type string "Klystron".
	 *	@return Type string "Klystron"
	 */
	virtual const std::string& GetType() const;

	/**
	 * Virtual constructor. Note a copy of a
	 * Klystron is attached to the same
	 * cavities in the model.
	 */
	virtual Klystron* Copy() const;

protected:

	std::vector<RFStructure*> rf_cavs;
	Mode kmode;
	double V_k;
	double phi_k;

	void AdjustCavities();
	void AdjustVectorSum(double v, double phi);
	void AppendBeamlineIndexes(std::vector<size_t>&) const;
};

/**
 * Class Klystron
 */
inline Klystron::~Klystron()
{
	// nothing to do
}

inline double Klystron::GetVoltage() const
{
	return V_k;
}

inline double Klystron::GetPhase() const
{
	return phi_k;
}

inline Complex Klystron::GetVoltagePhasor() const
{
	return Complex(V_k * cos(phi_k), V_k * sin(phi_k));
}

#endif
