/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _APERTURECONFIGURATION_H_
#define _APERTURECONFIGURATION_H_

#include <vector>
#include <string>
#include <fstream>

#include "AcceleratorModel.h"
#include "DataTable.h"

using namespace std;

class ApertureConfiguration
{
public:

	/**
	 * Constructor
	 */
	ApertureConfiguration();

	/**
	 * Constructor with an input file to load
	 * @param[in] InputFileName The name of the aperture file to load
	 */
	ApertureConfiguration(std::string InputFileName);

	/**
	 * Prints the configured apertures from AcceleratorModel
	 * @param[in] os The name of the stream to output to
	 */
	void OutputConfiguredAperture(AcceleratorModel* model, ostream& os);

	/**
	 * Configures the beam pipe for a given accelerator model
	 * @param[in] Model A pointer to the AcceleratorModel class to add the apertures to
	 */
	void ConfigureElementApertures(AcceleratorModel*);

	/**
	 * Deletes all apertures currently attached to the given accelerator model
	 * @param[in] Model A pointer to the AcceleratorModel class to add the apertures to
	 */
	void DeleteAllApertures(AcceleratorModel* Model);

	/**
	 * Set the stream for the collimator settings log.
	 * @param [in] os The stream to log the generated aperture to
	 */
	void SetLogFile(ostream& os);

	/**
	 * Enable/disable logging
	 * @param [in] flag The requested logging state
	 */
	void EnableLogging(bool flag);

	/**
	 * The output log file
	 */
	std::ostream* log;

	/**
	 * Enable/disable logging
	 */
	bool logFlag;

	/**
	 * The global list of Aperture entries
	 */
	DataTable ApertureDataTable;

};

#endif
