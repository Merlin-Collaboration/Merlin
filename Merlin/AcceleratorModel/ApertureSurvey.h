#ifndef ApertureSurvey_h
#define ApertureSurvey_h 1

#include <fstream>

#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/Aperture.h"

class ApertureSurvey
{
public:

	// Constructor
	ApertureSurvey(AcceleratorModel* model, std::string directory, double StepSize=0.1, size_t PointsPerElement=0);

	// Overloaded Constructor
	ApertureSurvey(AcceleratorModel* model, std::string directory, bool exact_s, double StepSize=0.1);

	// Overloaded Constructor
	ApertureSurvey(AcceleratorModel* model, std::ostream* os, double StepSize=0.1, size_t PointsPerElement=0);

	// Overloaded Constructor
	ApertureSurvey(AcceleratorModel* model, std::ostream* os, bool exact_s, double StepSize=0.1);

private:
	void CheckAperture(Aperture* ap, double s, double *lims);

	size_t Points;
	double StepSize;
	std::ofstream* OutputFile;
	AcceleratorModel* AccMod;
};

#endif

