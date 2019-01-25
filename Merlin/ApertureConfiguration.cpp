/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "AcceleratorModel.h"
#include "Collimator.h"
#include "DataTable.h"
#include "DataTableTFS.h"

#include "ApertureConfiguration.h"
#include "InterpolatedApertures.h"

using namespace std;

ApertureConfiguration::ApertureConfiguration() :
	log(nullptr), logFlag(false), ApertureEntry(nullptr), DefaultAperture(nullptr), DefaultApertureFlag(false)
{

}

ApertureConfiguration::ApertureConfiguration(string InputFileName) :
	log(nullptr), logFlag(false), ApertureEntry(nullptr), DefaultAperture(nullptr), DefaultApertureFlag(false)
{
	unique_ptr<DataTable> dt(DataTableReaderTFS(InputFileName).Read());
	unique_ptr<DataTable>& dtref = dt;
	AssignAperturesToList(dtref);
}

void ApertureConfiguration::AssignAperturesToList(unique_ptr<DataTable>& dt)
{
	ApertureFactory factory;
	for(auto row : *dt)
	{
		if(row.Get_d("APER_1") == 0 && row.Get_d("APER_2") == 0 && row.Get_d("APER_3") == 0
			&& row.Get_d("APER_4") == 0)
			continue;
		ApertureEntry = factory.getInstance(row.Get_s("APERTYPE"), row.Get_d("S") - row.Get_d("L"),
			row.Get_d("APER_1"), row.Get_d("APER_2"), row.Get_d("APER_3"), row.Get_d("APER_4"));
		ApertureList.push_back(ApertureEntry);
		ApertureEntry->setSlongitudinal(row.Get_d("S"));
		ApertureList.push_back(ApertureEntry);
	}
}

void ApertureConfiguration::OutputApertureList(std::ostream& os)
{
	cout << "APERTYPE" << "\t" << "S" << "\t" << "APER_1" << "\t" << "APER_2" << "\t" << "APER_3" << "\t" << "APER_4"
		 << endl;
	for(auto aperture : ApertureList)
	{
		cout << aperture->getApertureType() << "\t" << aperture->getSlongitudinal() << "\t"
			 << aperture->getRectHalfWidth() << "\t" << aperture->getRectHalfHeight() << "\t"
			 << aperture->getEllipHalfWidth() << "\t" << aperture->getEllipHalfHeight() << endl;
	}
}
void ApertureConfiguration::ConfigureElementApertures(AcceleratorModel* Model)
{
	std::vector<AcceleratorComponent*> Elements;
	int nElements = Model->ExtractTypedElements(Elements, "*");
	cout << "Got " << nElements << " elements for aperture configuration" << endl;
	cout << "Got " << ApertureList.size() << " Aperture entries" << endl;

	Aperture* aperture;
	ApertureFactory* factory;
	InterpolatorFactory* intfactory;

	for(std::vector<AcceleratorComponent*>::iterator comp = Elements.begin(); comp != Elements.end(); comp++)
	{
		if((*comp)->GetAperture() != nullptr)
			continue;
		if((*comp)->GetLength() == 0)
			continue;
		register double ElementLength = (*comp)->GetLength();
		register double Position = (*comp)->GetComponentLatticePosition();
		for(std::vector<Aperture*>::iterator itr = ApertureList.begin(); itr != ApertureList.end(); itr++)
		{
			if((*comp)->GetType() != "Drift" && (*itr)->getSlongitudinal() >= Position)
			{
				aperture = factory->getInstance((*itr)->getApertureType(), (*itr)->getSlongitudinal(),
					(*itr)->getRectHalfWidth(), (*itr)->getRectHalfHeight(), (*itr)->getEllipHalfWidth(),
					(*itr)->getEllipHalfHeight());
				(*comp)->SetAperture(aperture);
				itr = ApertureList.end();
				break;
			}
			if((*comp)->GetType() == "Drift"
				&& ((*itr)->getSlongitudinal() >= Position || itr == (ApertureList.end() - 1)))
			{
				std::vector<Aperture*> ThisElementAperture;
				if(itr == ApertureList.begin())
				{
					std::cout << "At first element " << (*comp)->GetQualifiedName()
							  << " getting aperture iterpolation from last element" << std::endl;

					//got the initial point
					itr = ApertureList.end() - 1;
					(*itr)->setSlongitudinal(0);
					ThisElementAperture.push_back(*itr);
					itr = ApertureList.begin();
				}
				else
				{
					itr--;
					ThisElementAperture.push_back(*itr);
					itr++;
				}
				while((*itr)->getSlongitudinal() <= (Position + ElementLength))
				{
					ThisElementAperture.push_back(*itr);
					itr++;
					if(itr == ApertureList.end())
					{
						itr = ApertureList.begin();
						(*itr)->setSlongitudinal(Position + ElementLength);
						ThisElementAperture.push_back(*itr);
						itr = ApertureList.end();
						break;
					}
				}
				if(itr != ApertureList.end())
				{
					ThisElementAperture.push_back(*itr);
					itr = ApertureList.end();
				}
				size_t NegativeCount = 0;
				for(size_t n = 0; n < ThisElementAperture.size(); n++)
				{
					if(ThisElementAperture[n]->getSlongitudinal() - Position < 0)
					{
						NegativeCount++;
					}
				}
				while(NegativeCount > 0)
				{
					ThisElementAperture.erase(ThisElementAperture.begin());
					NegativeCount--;
				}
				if(fequal(ThisElementAperture[0]->getSlongitudinal() - Position, 0.0, 1e-9)) //1e-7
				{
					ThisElementAperture[0]->setSlongitudinal(Position);
				}
				bool interpolate = false;
				bool typeChange = false;
				string type;

				for(size_t n = 1; n < ThisElementAperture.size(); n++)
				{
					if(ThisElementAperture[0]->getRectHalfWidth() != ThisElementAperture[n]->getRectHalfWidth()
						|| ThisElementAperture[0]->getRectHalfHeight()
						!= ThisElementAperture[n]->getRectHalfHeight()
						|| ThisElementAperture[0]->getEllipHalfWidth()
						!= ThisElementAperture[n]->getEllipHalfWidth()
						|| ThisElementAperture[0]->getEllipHalfHeight()
						!= ThisElementAperture[n]->getEllipHalfHeight())
					{
						interpolate = true;
					}
					else if(ThisElementAperture[0]->getType() != ThisElementAperture[n]->getType())
					{
						typeChange = true;
						if(ThisElementAperture[n]->getType() == "OCTAGON")
						{
							std::cerr << "Sorry, cannot interpolate from other geometries to an octagon" << std::endl;
							exit(EXIT_FAILURE);
						}
					}
				}
				if(interpolate)
				{
					bool octagon = false;
					vector<Aperture*> apVec;
					for(size_t n = 0; n < ThisElementAperture.size(); n++)
					{
						Aperture* apInterpolated = factory->getInstance(ThisElementAperture[n]->getType(),
							(ThisElementAperture[n]->getSlongitudinal() - Position),
							ThisElementAperture[n]->getRectHalfWidth(), ThisElementAperture[n]->getRectHalfHeight(),
							ThisElementAperture[n]->getEllipHalfWidth(),
							ThisElementAperture[n]->getEllipHalfHeight());
						apVec.push_back(apInterpolated);
					}
					aperture = intfactory->getInstance(apVec);
					(*comp)->SetAperture(aperture);
					itr = ApertureList.end();
					break;
				}
			}
			else
			{
				aperture = factory->getInstance((*itr)->getApertureType(), (*itr)->getSlongitudinal(),
					(*itr)->getRectHalfWidth(), (*itr)->getRectHalfHeight(), (*itr)->getEllipHalfWidth(),
					(*itr)->getEllipHalfHeight());
				(*comp)->SetAperture(aperture);
				itr = ApertureList.end();
				break;
			}
		}
		if(logFlag)
		{
			*log << setw(25) << left << (*comp)->GetName();
			*log << setw(14) << left << (*comp)->GetType();
			*log << setw(10) << left << (*comp)->GetLength();
			*log << setw(10) << left << (*comp)->GetComponentLatticePosition();
			if((*comp)->GetAperture() != nullptr)
			{
				(*comp)->GetAperture()->printout(*log);
			}
			*log << endl;
		}
	}
}

void ApertureConfiguration::SetLogFile(ostream& os)
{
	log = &os;
}

void ApertureConfiguration::EnableLogging(bool flg)
{
	logFlag = flg;
}

void ApertureConfiguration::DeleteAllApertures(AcceleratorModel* Model)
{
	std::vector<AcceleratorComponent*> Elements;
	int nElements = Model->ExtractTypedElements(Elements, "*");

	for(std::vector<AcceleratorComponent*>::iterator comp = Elements.begin(); comp != Elements.end(); comp++)
	{
		if((*comp)->GetAperture() != nullptr)
		{
			delete (*comp)->GetAperture();
			(*comp)->SetAperture(nullptr);
		}
	}
}

void ApertureConfiguration::SetDefaultAperture(Aperture* ap)
{
	DefaultAperture = ap;
}

void ApertureConfiguration::EnableDefaultAperture(bool flag)
{
	DefaultApertureFlag = flag;
}
