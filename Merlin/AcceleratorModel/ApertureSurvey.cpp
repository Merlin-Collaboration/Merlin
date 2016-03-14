#include <iostream>
#include <sstream>
#include <string>

#include "AcceleratorModel/ApertureSurvey.h"
#include "AcceleratorModel/StdComponent/Collimator.h"

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"


ApertureSurvey::ApertureSurvey(AcceleratorModel* Model, std::string Directory, double Step, size_t PointsPerElement)
{
	StepSize = Step;
	Points = PointsPerElement;
	AccMod = Model;

	//Create a file in the given directory
	std::ostringstream FileStreamName;
	FileStreamName << Directory << "Aperture_Survey_" << Step << "_steps_OR_" << Points << "_points.txt";
	std::ofstream* OutputFile = new std::ofstream(FileStreamName.str().c_str());

	if(!OutputFile->good())
	{
		std::cerr << "Could not open ApertureSurvey file" << std::endl;
		exit(EXIT_FAILURE);
	}

	(*OutputFile) << "#name\ttype\ts_end\tlength\tap_px\tap_mx\tap_py\tap_my" << std::endl;

	double s = 0;
	double LastSample = 0-StepSize;
	double lims[4];

	for (AcceleratorModel::BeamlineIterator bi = AccMod->GetBeamline().begin(); bi != AccMod->GetBeamline().end(); bi++)
	{
		AcceleratorComponent *ac = &(*bi)->GetComponent();

		if (fabs(s - ac->GetComponentLatticePosition())> 1e-6)
		{
			exit(EXIT_FAILURE);
		}

		Aperture* ap =	ac->GetAperture();
		Collimator* aCollimator = dynamic_cast<Collimator*>(ac);

		std::vector<double> zs;
		if (Points > 0)
		{
			for (size_t i=0; i<Points; i++)
			{
				zs.push_back(i*ac->GetLength()*(1.0/(Points-1)));
			}
		}
		else
		{
			while (LastSample+StepSize < s + ac->GetLength())
			{
				LastSample += StepSize;
				zs.push_back(LastSample-s);
			}
		}

		for (size_t zi = 0; zi < zs.size(); zi++)
		{
			double z = zs[zi];

			if (ap != nullptr)
			{
				ApertureSurvey::CheckAperture(ap, z, lims);
			}
			else
			{
				lims[0] = lims[1] = lims[2]= lims[3] = 1;
			}

			(*OutputFile) << ac->GetName() << "\t";
			(*OutputFile) << ac->GetType() << "\t";
			(*OutputFile) << ac->GetComponentLatticePosition()+ac->GetLength() << "\t";
			(*OutputFile) << ac->GetLength() << "\t";
			//~ (*os) << s+z << "\t";
			(*OutputFile) << lims[0] << "\t";
			(*OutputFile) << lims[1] << "\t";
			(*OutputFile) << lims[2] << "\t";
			(*OutputFile) << lims[3] << std::endl;
		}
		s += ac->GetLength();
	}
}

ApertureSurvey::ApertureSurvey(AcceleratorModel* Model, std::string Directory, bool exact_s, double Step)
{
	if(!exact_s)
	{
		ApertureSurvey(Model, Directory, Step);
	}
	else
	{
		StepSize = Step;
		Points = 0;
		AccMod = Model;

		//Create a file in the given directory
		std::ostringstream FileStreamName;
		FileStreamName << Directory << "Aperture_Survey_" << Step << "_steps_OR_" << Points << "_points.txt";
		std::ofstream* OutputFile = new std::ofstream(FileStreamName.str().c_str());

		if(!OutputFile->good())
		{
			std::cerr << "Could not open ApertureSurvey file" << std::endl;
			exit(EXIT_FAILURE);
		}

		(*OutputFile) << "#name\ttype\ts_end\tlength\tap_px\tap_mx\tap_py\tap_my" << std::endl;

		double s = 0;
		double LastSample = 0-StepSize;
		double lims[4];

		for (AcceleratorModel::BeamlineIterator bi = AccMod->GetBeamline().begin(); bi != AccMod->GetBeamline().end(); bi++)
		{
			AcceleratorComponent *ac = &(*bi)->GetComponent();

			if (fabs(s - ac->GetComponentLatticePosition()) > 1e-6)
			{
				exit(EXIT_FAILURE);
			}

			Aperture* ap =	ac->GetAperture();
			Collimator* aCollimator = dynamic_cast<Collimator*>(ac);

			std::vector<double> zs;
			if (Points > 0)
			{
				for (size_t i=0; i<Points; i++)
				{
					zs.push_back(i*ac->GetLength()*(1.0/(Points-1)));
				}
			}
			else
			{
				while (LastSample+StepSize < s + ac->GetLength())
				{
					LastSample += StepSize;

					zs.push_back(LastSample-s);
				}
			}
			for (size_t zi = 0; zi < zs.size(); zi++)
			{
				double z = zs[zi];

				if (ap != nullptr)
				{
					ApertureSurvey::CheckAperture(ap, z, lims);
				}
				else
				{
					lims[0] = lims[1] = lims[2]= lims[3] = 1;
				}

				(*OutputFile) << ac->GetName() << "\t";
				(*OutputFile) << ac->GetType() << "\t";
				(*OutputFile) << ac->GetComponentLatticePosition()+z << "\t";
				//~ (*OutputFile) << ac->GetComponentLatticePosition()+ac->GetLength() << "\t";
				(*OutputFile) << ac->GetLength() << "\t";
				(*OutputFile) << lims[0] << "\t";
				(*OutputFile) << lims[1] << "\t";
				(*OutputFile) << lims[2] << "\t";
				(*OutputFile) << lims[3] << std::endl;
			}
			s += ac->GetLength();
		}
	}
}

ApertureSurvey::ApertureSurvey(AcceleratorModel* Model, std::ostream* os, double Step, size_t PointsPerElement)
{
	StepSize = Step;
	Points = PointsPerElement;
	AccMod = Model;

	(*os) << "#name\ttype\ts_end\tlength\tap_px\tap_mx\tap_py\tap_my" << std::endl;

	double s = 0;
	double LastSample = 0-StepSize;
	double lims[4];

	for (AcceleratorModel::BeamlineIterator bi = AccMod->GetBeamline().begin(); bi != AccMod->GetBeamline().end(); bi++)
	{
		AcceleratorComponent *ac = &(*bi)->GetComponent();

		if (fabs(s - ac->GetComponentLatticePosition())> 1e-6)
		{
			exit(EXIT_FAILURE);
		}

		Aperture* ap =	ac->GetAperture();
		Collimator* aCollimator = dynamic_cast<Collimator*>(ac);

		vector<double> zs;
		if (Points > 0)
		{
			for (size_t i=0; i<Points; i++)
			{
				zs.push_back(i*ac->GetLength()*(1.0/(Points-1)));
			}
		}
		else
		{
			while (LastSample+StepSize < s + ac->GetLength())
			{
				LastSample += StepSize;

				zs.push_back(LastSample-s);
			}
		}
		for (size_t zi = 0; zi < zs.size(); zi++)
		{
			double z = zs[zi];

			if (ap != nullptr)
			{
				ApertureSurvey::CheckAperture(ap, z, lims);
			}
			else
			{
				lims[0] = lims[1] = lims[2] = lims[3] = 1;
			}

			(*os) << ac->GetName() << "\t";
			(*os) << ac->GetType() << "\t";
			(*os) << ac->GetComponentLatticePosition()+ac->GetLength() << "\t";
			(*os) << ac->GetLength() << "\t";
			//~ (*os) << s+z << "\t";
			(*os) << lims[0] << "\t";
			(*os) << lims[1] << "\t";
			(*os) << lims[2] << "\t";
			(*os) << lims[3] << std::endl;
		}
		s += ac->GetLength();
	}
}

ApertureSurvey::ApertureSurvey(AcceleratorModel* Model, std::ostream* os, bool exact_s, double Step)
{
	if(!exact_s)
	{
		ApertureSurvey(Model, os, Step);
	}
	else
	{
		StepSize = Step;
		Points = 0;
		AccMod = Model;

		(*os) << "#name\ttype\ts_end\tlength\tap_px\tap_mx\tap_py\tap_my" << std::endl;

		double s = 0;
		double LastSample = 0-StepSize;
		double lims[4];

		for (AcceleratorModel::BeamlineIterator bi = AccMod->GetBeamline().begin(); bi != AccMod->GetBeamline().end(); bi++)
		{

			AcceleratorComponent *ac = &(*bi)->GetComponent();

			if (fabs(s - ac->GetComponentLatticePosition())> 1e-6)
			{
				exit(EXIT_FAILURE);
			}

			Aperture* ap =	ac->GetAperture();
			Collimator* aCollimator = dynamic_cast<Collimator*>(ac);

			std::vector<double> zs;

			if (Points > 0)
			{
				for (size_t i=0; i<Points; i++)
				{
					zs.push_back(i*ac->GetLength()*(1.0/(Points-1)));
				}
			}
			else
			{
				while (LastSample+StepSize < s + ac->GetLength())
				{
					LastSample += StepSize;
					zs.push_back(LastSample-s);
				}
			}

			for (size_t zi = 0; zi < zs.size(); zi++)
			{
				double z = zs[zi];

				if (ap != nullptr)
				{
					ApertureSurvey::CheckAperture(ap, z, lims);
				}
				else
				{
					lims[0] = lims[1] = lims[2] = lims[3] = 1;
				}

				(*os) << ac->GetName() << "\t";
				(*os) << ac->GetType() << "\t";
				(*os) << ac->GetComponentLatticePosition()+z << "\t";
				//~ (*os) << ac->GetComponentLatticePosition()+ac->GetLength() << "\t";
				(*os) << ac->GetLength() << "\t";
				(*os) << lims[0] << "\t";
				(*os) << lims[1] << "\t";
				(*os) << lims[2] << "\t";
				(*os) << lims[3] << std::endl;
			}
			s += ac->GetLength();
		}
	}
}

void ApertureSurvey::CheckAperture(Aperture* ap, double s, double *aps)
{
	const double Step = 1e-6;
	const double max = 1.0;
	const double min = 0.0;

	// iterate through directions
	for(int dir=0; dir<4; dir++)
	{
		double xdir=0, ydir=0;

		if(dir==0)
		{
			xdir=+1;
		}
		else if (dir==1)
		{
			xdir=-1;
		}
		else if (dir==2)
		{
			ydir=+1;
		}
		else if (dir==3)
		{
			ydir=-1;
		}

		aps[dir] = 0;

		// scan for limit
		double below=min, above=max;

		while(above-below > Step)
		{
			double guess = (above+below)/2;

			if (ap->PointInside(xdir*guess, ydir*guess, s))
			{
				below = guess;
			}
			else
			{
				above = guess;
			}
		}
		aps[dir] = (above+below)/2;
	}
}

