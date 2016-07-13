#ifndef ApertureSurvey_h
#define ApertureSurvey_h 1

#include <string>
#include <iostream>

#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/Aperture.h"

namespace ApertureSurvey
{

// Constructor
void ApertureSurvey(AcceleratorModel* model, string file_name, double step_size=0.1, size_t points_per_element=0);
// Overloaded Constructor
void ApertureSurvey(AcceleratorModel* model, string file_name, bool exact_s, double step_size=0.1);
// Overloaded Constructor
void ApertureSurvey(AcceleratorModel* model, std::ostream* os, double step_size=0.1, size_t points_per_element=0);
// Overloaded Constructor
void ApertureSurvey(AcceleratorModel* model, std::ostream* os, bool exact_s, double step_size=0.1);

void CheckAperture(Aperture* ap, double s, double *aps);

}

#endif
