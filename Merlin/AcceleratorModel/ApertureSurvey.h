#ifndef ApertureSurvey_h
#define ApertureSurvey_h 1

#include <fstream>

#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/Aperture.h"

class ApertureSurvey
{
public:
	// Constructor
	ApertureSurvey(AcceleratorModel* model, string directory, double step_size=0.1, size_t points_per_element=0);
	// Overloaded Constructor
	ApertureSurvey(AcceleratorModel* model, string directory, bool exact_s, double step_size=0.1);
	// Overloaded Constructor
	ApertureSurvey(AcceleratorModel* model, std::ostream* os, double step_size=0.1, size_t points_per_element=0);
	// Overloaded Constructor
	ApertureSurvey(AcceleratorModel* model, std::ostream* os, bool exact_s, double step_size=0.1);
	
private:
	void CheckAperture(Aperture* ap, double s, double *lims);
	
	size_t points;
	double step_size;	
	std::ofstream* output_file;
	AcceleratorModel* AccMod;
};

#endif
