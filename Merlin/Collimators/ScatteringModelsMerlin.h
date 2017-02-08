#ifndef ScatteringModelsMerlin_h
#define ScatteringModelsMerlin_h 1

#include "ScatteringModel.h"

namespace Collimation
{
/**
 * Base class for fixed scattering models
 *
 * Classes derived from this should set up their scattering processes in
 * their constructor, and then set is_fixed to true.
 */
class ScatteringModelFixed: public ScatteringModel
{
public:
	ScatteringModelFixed();
	virtual ~ScatteringModelFixed();

	/**
	 * For fixed scattering models, this cannot be called from outside
	 * of the constructor.
	 */
	void AddProcess(Collimation::ScatteringProcess* S);

protected:
	bool is_fixed;
};

/**
 * Merlin physics model
 *
 * ScatteringModel preset with Merlin physics model as described in
 * Appleby, R.B., Barlow, R.J., Molson, J.G. et al.
 * Eur. Phys. J. C (2016) 76: 520.
 * https://dx.doi.org/10.1140/epjc/s10052-016-4363-7
 */
class ScatteringModelMerlin: public ScatteringModelFixed
{
public:
	ScatteringModelMerlin();
};

}
#endif
