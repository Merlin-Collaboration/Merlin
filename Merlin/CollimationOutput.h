/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CollimationOutput_h
#define CollimationOutput_h 1

#include <string>
#include <vector>

#include "AcceleratorComponent.h"

#include "ParticleBunch.h"

#include "PSTypes.h"

namespace ParticleTracking
{

/**
 * Struct used to store individual lost particle data
 */
struct LossData
{
	typedef enum
	{
		Collimator,
		Cold,
		Warm,
		Undefined

	} LossTypes;

	std::string ElementName;
	PSvector p;
	double s;
	double interval;
	double position;
	double length;
	double lost;
	LossTypes temperature;
	int turn;
	int coll_id;
	double angle;

	LossData() :
		ElementName(), p(), s(), interval(), position(), length(), lost(), turn(), coll_id(), angle()
	{
	}

	void reset()
	{
		ElementName = "_";
		s = 0;
		interval = 0;
		position = 0;
		length = 0;
		lost = 0;
		temperature = Undefined;
		turn = 0;
		coll_id = 0;
		angle = 0;
	}

	bool operator<(LossData other) const
	{
		return (s + position) > (other.s + other.position);
	}

	bool operator==(LossData other) const
	{
		if((s + position) == (other.s + other.position))
		{
			return true;
		}

		return false;
	}

	/**
	 * Note that the + operator cannot preserve the particle PSvector p
	 */
	LossData operator+(LossData other)
	{
		/**
		 * Create temporary LossData struct to hold final LossData object
		 */
		LossData temp;

		/**
		 * Check that the loss is in the same element
		 */
		if(ElementName == other.ElementName)
		{
			temp.ElementName = ElementName;
			temp.p = p;
			temp.s = s;
			temp.interval = interval;
			temp.position = position;
			temp.length = length;
			temp.temperature = temperature;
			temp.lost = lost + other.lost;

			return temp;
		}
		else
		{
			std::cout
				<<
				"Warning: CollimationOutput Class: Cannot operator+ for losses in different elements, returning original LossData object"
				<< std::endl;
			return *this;
		}
	}

	LossData operator++()
	{
		lost += 1;
		return *this;
	}

};

// Comparison function used to sort losses in order of s position
inline bool Compare_LossData(const LossData &a, const LossData &b)
{
	return (a.s + a.position + a.interval) < (b.s + b.position + a.interval);
}

inline bool Merge_LossData(const LossData &a, const LossData &b)
{
	if((a.s + a.position + a.interval) == (b.s + b.position + a.interval))
	{
		return true;
	}
	return false;
}

// Possible output types for each class
typedef enum
{
	nearestelement,
	precise,
	tencm

} OutputType;

/**
 * CollimationOutput handles the output from the collimation process, specifically
 * lost particles. It is called from CollimateProtonProcess::DoScatter and
 * allows the user to create loss map output files, root hist files, or
 * a user specified output format.
 */
class CollimationOutput
{

public:

	/**
	 * Constructor
	 */
	CollimationOutput(OutputType otype = nearestelement);

	/**
	 * Destructor
	 */
	~CollimationOutput();

	/**
	 * Finalise will call any sorting algorithms and perform formatting for final output
	 */
	virtual void Finalise()
	{
	}

	/**
	 * Perform the final output
	 */
	virtual void Output(std::ostream* os)
	{
	}

	/**
	 * Called from CollimateProtonProcess::DoScatter to add a particle to the CollimationOutput
	 */
	virtual void Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn = 0)
	{
	}

	/**
	 * Output type switch
	 */
	OutputType otype;

	/**
	 * Temporary LossData struct use to transfer data
	 */
	LossData temp;

	/**
	 * Vector to hold the loss data
	 */
	std::vector<LossData> DeadParticles;

	/**
	 * Vector to hold output data
	 */
	std::vector<LossData> OutputLosses;

protected:
	AcceleratorComponent* currentComponent;

private:
};

} //End namespace ParticleTracking

#endif
