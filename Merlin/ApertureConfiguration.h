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
	 * Load the Aperture settings from an input file.
	 * @param[in] InputFileName The name of the aperture file to load
	 */
	void AssignAperturesToList(unique_ptr<DataTable>&);

	/**
	 * Dumps the input file back out
	 * @param[in] os The name of the stream to output to
	 */
	void OutputApertureList(std::ostream& os);

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
	 * Set a default class of aperture to use in ambiguous cases
	 * @param [in] flag The requested logging state
	 */
	void SetDefaultAperture(Aperture* ap);

	/**
	 * Enable/disable use of the default aperture where it is not possible to clearly select an aperture type (e.g. OCTAGON -> RECTELLIPSE joins)
	 * @param [in] flag A bool to enable or disable the usage of the default aperture
	 */
	void EnableDefaultAperture(bool flag);

	/**
	 * The output log file
	 */
	std::ostream* log;

	/**
	 * Enable/disable logging
	 */
	bool logFlag;

	/**
	 * Typedef for access to the enum
	 */
	typedef size_t ApertureClass_t;

	/**
	 * One aperture entry
	 */
	Aperture* ApertureEntry;
	//ap ApertureEntry;
	/**
	 * The global list of Aperture entries
	 */
	vector<Aperture*> ApertureList;
	//std::vector<ap> ApertureList;

	/**
	 * A pointer to a default aperture entry
	 */
	Aperture* DefaultAperture;
	bool DefaultApertureFlag;
};

#endif
