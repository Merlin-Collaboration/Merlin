/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ApertureSurvey_h
#define ApertureSurvey_h 1

#include <string>
#include <iostream>

#include "AcceleratorModel.h"
#include "Aperture.h"

/**
 * Survey an accelerator and find the aperture at each step.
 *
 * At each point the aperture is measured by bisection search, using the
 * PointInside() method of the aperture attached to the element. This
 * makes it useful for debugging the whole aperture system.
 *
 * Output file contains the columns:
 *
 * The element's
 * name, type, s_end, length
 *
 * Position
 * s
 *
 * Aperture limit in plus and minus x and y directions
 * ap_px ap_mx ap_py ap_my
 */
namespace ApertureSurvey
{
//! Mode for survey. Changes how step is interpreted.
enum SurveyType
{
	points_per_element,
	/**< Record \a step times per element, e.g. 1: just at start, 2: start and end, 3: start, end and middle. */
	distance,
	/**< Record every \a step distance, within each element. Include start and end of each element. */
	abs_distance /**< Record every \a step distance, along \a model. */

};

/**
 * Run aperture survey
 * @param model The accelerator mode
 * @param os Output stream
 * @param mode Stepping mode
 * @param step Definition depends on \a mode
 */
void ApertureSurvey(AcceleratorModel* model, std::ostream* os, SurveyType mode, double step);

/**
 * Run aperture survey
 * @param model The accelerator mode
 * @param filename Filename to output to
 * @param mode Stepping mode
 * @param step Definition depends on \a mode
 */
void ApertureSurvey(AcceleratorModel* model, string file_name, SurveyType mode, double step);

}

#endif
