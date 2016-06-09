#include "Collimators/Output/FlukaCollimationOutput.h"

#include "AcceleratorModel/AcceleratorComponent.h"
#include "AcceleratorModel/StdComponent/Collimator.h"
#include "AcceleratorModel/Apertures/CollimatorAperture.h"

namespace ParticleTracking
{

FlukaCollimationOutput::FlukaCollimationOutput(OutputType ot)
{
	otype = ot;
}

void FlukaCollimationOutput::Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn)
{
	// If current component is a collimator we store the loss, otherwise we do not
	if (currentComponent != &currcomponent)
	{
		currentComponent = &currcomponent;
	}
	temp.reset();

	Collimator* aCollimator = dynamic_cast<Collimator*>(&currcomponent);
	bool is_collimator = aCollimator;

	if(is_collimator)
	{

		temp.ElementName = currentComponent->GetQualifiedName().c_str();
		temp.s = currentComponent->GetComponentLatticePosition();

		// For Fluka output pos is the lost position within the element
		temp.position = pos;
		temp.length = currentComponent->GetLength();
		temp.lost = 1;

		temp.coll_id = currentComponent->GetCollID();
		temp.turn = turn;

		const CollimatorAperture* tap= dynamic_cast<const CollimatorAperture*> (currentComponent->GetAperture());
		temp.angle = tap->GetCollimatorTilt();

		temp.p = particle;

		//pushback vector
		DeadParticles.push_back(temp);
	}

}

void FlukaCollimationOutput::Finalise()
{
	for(std::vector <LossData>::iterator its = DeadParticles.begin(); its != DeadParticles.end(); ++its)
	{
		if( (*its).p.type() == 1 || (*its).p.type() == 4 )
		{
			OutputLosses.push_back(*its);
		}
	}
}

void FlukaCollimationOutput::Output(std::ostream* os)
{
	cout << "\nFlukaCollimationOutput OutputLosses size = " << OutputLosses.size() << ", DeadParticles.size() = " << DeadParticles.size() << endl;
	(*os) << "#\t1=icoll\t2=c_rotation\t3=s\t4=x\t5=xp\t6=y\t7=yp\t8=nabs\t9=np\t10=ntu" << endl;
	for(std::vector <LossData>::iterator its = OutputLosses.begin(); its != OutputLosses.end(); ++its)
	{
		(*os) << setw(16) << left << (*its).coll_id;
		(*os) << setw(20) << left << (*its).angle;
		(*os) << setw(20) << left << (*its).position;
		(*os) << setw(20) << left << (*its).p.x();
		(*os) << setw(20) << left << (*its).p.xp();
		(*os) << setw(20) << left << (*its).p.y();
		(*os) << setw(20) << left << (*its).p.yp();
		(*os) << setw(20) << left << (*its).p.type();
		(*os) << setw(20) << left << (*its).p.id();
		(*os) << setw(20) << left << (*its).turn;
		(*os) << endl;
	}
}

}//End namespace

