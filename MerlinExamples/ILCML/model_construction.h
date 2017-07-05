#ifndef _h_model_construction
#define _h_model_construction

#include "merlin_config.h"
#include "AcceleratorModel.h"
#include "ComponentFrame.h"
#include "BeamData.h"
#include <utility>

std::pair<AcceleratorModel*,BeamData*> ConstructModel(const string& fname);

#endif



