/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <iomanip>

#include "AcceleratorModel.h"
#include "Collimator.h"
#include "DataTable.h"
#include "DataTableTFS.h"

#include "ApertureConfiguration.h"
#include "InterpolatedApertures.h"

using namespace std;

ApertureConfiguration::ApertureConfiguration(string InputFileName)
{
	unique_ptr<DataTable> dt;
	try
	{
		dt = make_unique<DataTable>(DataTableReaderTFS(InputFileName).Read());
	}
	catch(BadFormatException& e)
	{
		std::cerr << "ApertureConfiguration: Error reading " << InputFileName << std::endl;
		throw;
	}

	ApertureDataTable.AddColumn("S", 'd');
	ApertureDataTable.AddColumn("APER_1", 'd');
	ApertureDataTable.AddColumn("APER_2", 'd');
	ApertureDataTable.AddColumn("APER_3", 'd');
	ApertureDataTable.AddColumn("APER_4", 'd');
	ApertureDataTable.AddColumn("APERTYPE", 's');
	//filter non-aperture entries
	size_t row = 0;
	for(auto &itr : *dt)
	{
		if(itr.Get_d("APER_1") != 0 || itr.Get_d("APER_2") != 0 || itr.Get_d("APER_3") != 0 || itr.Get_d(
				"APER_4") != 0)
		{
			ApertureDataTable.ApertureFromRow(itr);
			ApertureDataTable.Set("S", row, (itr.Get_d("S") - itr.Get_d("L")));
			++row;
			ApertureDataTable.ApertureFromRow(itr);
			++row;
		}
	}
}

void ApertureConfiguration::ConfigureElementApertures(AcceleratorModel* Model)
{
	std::vector<AcceleratorComponent*> Elements;
	int nElements = Model->ExtractTypedElements(Elements, "*");
	cout << "Got " << nElements << " elements for aperture configuration" << endl;
	cout << "Got " << ApertureDataTable.Length() << " Aperture entries" << endl;
	size_t noConfiguredApertures = 0;

	Aperture* aperture;

	for(auto &comp :Elements)
	{
		if(comp->GetAperture() != nullptr)
			continue;
		if(comp->GetLength() == 0)
			continue;
		double ElementLength = comp->GetLength();
		double Position = comp->GetComponentLatticePosition();
		if(comp->GetType() != "Drift")
		{
			auto this_ap = std::lower_bound(ApertureDataTable.begin(), ApertureDataTable.end(), Position,
				[](const auto a, const auto b){
				return a.Get_d("S") < b;
			});

			aperture = ApertureFactory::GetInstance(*this_ap);
			comp->SetAperture(aperture);
			++noConfiguredApertures;
		}
		else
		{
			auto this_ap = std::lower_bound(ApertureDataTable.begin(), ApertureDataTable.end(), Position,
				[](const auto a, const auto b){
				return a.Get_d("S") < b;
			});

			if(this_ap == ApertureDataTable.end())
			{
				--this_ap;
			}

			auto last_ap = ApertureDataTable.end() - 1;

			while(this_ap != ApertureDataTable.end())
			{
				size_t row = 0;
				DataTable ThisElementAperture;
				DataTable CleanElementAperture;
				ThisElementAperture.AddColumn("APERTYPE", 's');
				ThisElementAperture.AddColumn("S", 'd');
				ThisElementAperture.AddColumn("APER_1", 'd');
				ThisElementAperture.AddColumn("APER_2", 'd');
				ThisElementAperture.AddColumn("APER_3", 'd');
				ThisElementAperture.AddColumn("APER_4", 'd');

				CleanElementAperture = ThisElementAperture;

				if(this_ap == ApertureDataTable.begin())
				{
					//get initial point, interpolate from last point
					ThisElementAperture.ApertureFromRow(*last_ap);
					++row;
				}
				else
				{
					--this_ap;
					ThisElementAperture.ApertureFromRow(*this_ap);
					++this_ap;
				}
				while((*this_ap).Get_d("S") <= (Position + ElementLength))
				{
					++row;
					ThisElementAperture.ApertureFromRow(*this_ap);
					++this_ap;
					if(this_ap == ApertureDataTable.end())
					{
						break;
					}
				}
				if(this_ap == ApertureDataTable.end())
				{
					++row;
					this_ap = ApertureDataTable.begin();
					ThisElementAperture.ApertureFromRow(*this_ap);
				}
				else
				{
					++row;
					ThisElementAperture.ApertureFromRow(*this_ap);
					this_ap = ApertureDataTable.end();
				}

				//First do a little bit of cleaning
				//If we have an entry at 0 (or very close to), and also an entry at negative values, we can discard the negative entry
				size_t thisit = 0;
				size_t cleanit = 0;
				bool zeroentry = false;
				size_t negcount = 0;

				for(const auto itr : ThisElementAperture)
				{
					if(itr == *ThisElementAperture.begin())
					{
						for(size_t row = 0; row < ThisElementAperture.Length(); ++row)
						{
							if(ThisElementAperture.Get_d("S", row) - Position < 0)
							{
								++negcount;
							}
							if(fequal(ThisElementAperture.Get_d("S", row) - Position, 0.0, 1e-7))
							{
								ThisElementAperture.Set("S", row, Position);
								zeroentry = true;
							}
						}
					}
					if(negcount != 0 && zeroentry == true)
					{
						--negcount;
						continue;
					}
					if(negcount > 1)
					{
						--negcount;
						continue;
					}
					CleanElementAperture.ApertureFromRow(itr);
					++cleanit;
				}

				if(fequal(CleanElementAperture.Get_d("S", 0) - Position, 0.0, 5e-7))
				{
					CleanElementAperture.Set("S", 0, Position);
				}

				//Determine need to interpolate
				bool interpolate = false;
				bool typeChange = false;
				for(auto &itr :CleanElementAperture)
				{
					if(CleanElementAperture.begin()->Get_d("APER_1") != itr.Get_d("APER_1")
						|| CleanElementAperture.begin()->Get_d("APER_2") != itr.Get_d("APER_2")
						|| CleanElementAperture.begin()->Get_d("APER_3") != itr.Get_d("APER_3")
						|| CleanElementAperture.begin()->Get_d("APER_4") != itr.Get_d("APER_4"))
					{
						interpolate = true;
					}
					if(typeChange == false && CleanElementAperture.begin()->Get_s("APERTYPE") != itr.Get_s(
							"APERTYPE"))
					{
						typeChange = true;
					}
				}
				if(interpolate)
				{
					for(size_t row = 0; row < CleanElementAperture.Length(); ++row)
					{
						CleanElementAperture.Set_d("S", row, (CleanElementAperture.Get_d("S", row) - Position));
					}
					if(typeChange == false)
					{
						aperture = InterpolatorFactory::GetInstance(CleanElementAperture);
						comp->SetAperture(aperture);
						++noConfiguredApertures;
						break;
					}
					else
					{
						bool rectellipse = false;
						bool octagon = false;
						for(size_t row = 0; row < CleanElementAperture.Length(); ++row)
						{
							//CHANGE ZERO VALUE EXCEPTIONS
							if(CleanElementAperture.Get_s("APERTYPE", row) == "RECTELLIPSE")
							{
								CleanElementAperture.Set("APER_1", row, CleanElementAperture.Get_d("APER_1", row));
								CleanElementAperture.Set("APER_2", row, CleanElementAperture.Get_d("APER_2", row));
								CleanElementAperture.Set("APER_3", row, CleanElementAperture.Get_d("APER_3", row));
								CleanElementAperture.Set("APER_4", row, CleanElementAperture.Get_d("APER_4", row));
								rectellipse = true;
							}
							if(CleanElementAperture.Get_s("APERTYPE", row) == "CIRCLE")
							{
								CleanElementAperture.Set("APER_1", row, CleanElementAperture.Get_d("APER_1", row));
								CleanElementAperture.Set("APER_2", row, CleanElementAperture.Get_d("APER_1", row));
								CleanElementAperture.Set("APER_3", row, CleanElementAperture.Get_d("APER_1", row));
								CleanElementAperture.Set("APER_4", row, CleanElementAperture.Get_d("APER_1", row));
								rectellipse = true;
							}
							if(CleanElementAperture.Get_s("APERTYPE", row) == "ELLIPSE")
							{
								CleanElementAperture.Set("APER_1", row, CleanElementAperture.Get_d("APER_1", row));
								CleanElementAperture.Set("APER_2", row, CleanElementAperture.Get_d("APER_2", row));
								CleanElementAperture.Set("APER_3", row, CleanElementAperture.Get_d("APER_1", row));
								CleanElementAperture.Set("APER_4", row, CleanElementAperture.Get_d("APER_2", row));
								rectellipse = true;
							}
							if(CleanElementAperture.Get_s("APERTYPE", row) == "RECTANGLE")
							{
								CleanElementAperture.Set("APER_1", row, CleanElementAperture.Get_d("APER_1", row));
								CleanElementAperture.Set("APER_2", row, CleanElementAperture.Get_d("APER_2", row));
								CleanElementAperture.Set("APER_3", row, CleanElementAperture.Get_d("APER_1", row));
								CleanElementAperture.Set("APER_4", row, CleanElementAperture.Get_d("APER_2", row));
							}
							if(CleanElementAperture.Get_s("APERTYPE", row) == "OCTAGON")
							{
								CleanElementAperture.Set("APER_1", row, CleanElementAperture.Get_d("APER_1", row));
								CleanElementAperture.Set("APER_2", row, CleanElementAperture.Get_d("APER_2", row));
								CleanElementAperture.Set("APER_3", row, CleanElementAperture.Get_d("APER_3", row));
								CleanElementAperture.Set("APER_4", row, CleanElementAperture.Get_d("APER_4", row));
								octagon = true;
							}
						}
						if(rectellipse && octagon)
						{
							cerr
								<<
								"Attempting to interpolate incompatible aperture types (RECTELLIPSE and OCTAGON)... exiting simulation."
								<< endl;
							exit(EXIT_FAILURE);
						}
						aperture = InterpolatorFactory::GetInstance(CleanElementAperture);
						comp->SetAperture(aperture);
						++noConfiguredApertures;
						break;
					}
				}
				else
				{
					DataTableRowIterator itr = CleanElementAperture.begin();
					aperture = ApertureFactory::GetInstance(*itr);
					comp->SetAperture(aperture);
					++noConfiguredApertures;
					break;
				}
				this_ap++;
			}
		}
		if(log)
		{
			*log << setw(25) << left << comp->GetName();
			*log << setw(14) << left << comp->GetType();
			*log << setw(10) << left << comp->GetLength();
			*log << setw(10) << left << comp->GetComponentLatticePosition();
			if(comp->GetAperture() != nullptr)
			{
				comp->GetAperture()->printout(*log);
			}
			*log << endl;
		}
	}
	cout << "Number of configured apertures: " << noConfiguredApertures << endl;
}

void ApertureConfiguration::OutputConfiguredAperture(AcceleratorModel* Model, ostream& os)
{
	std::vector<AcceleratorComponent*> Elements;
	int nElements = Model->ExtractTypedElements(Elements, "*");
	size_t number = 1;

	for(const auto &comp : Elements)
	{
		if(comp->GetAperture() != nullptr)
		{
			cout << number << "\t"
				 << comp->GetIndex() << "\t"
				 << comp->GetName() << "\t"
				 << comp->GetQualifiedName() << "\t"
				 << comp->GetEMField() << "\t"
				 << comp->GetType() << "\t"
				 << comp->GetComponentLatticePosition() << "\t"
				 << comp->GetGeometry()->GetGeometryLength() << "\t"
				 << comp->GetIndex() << "\t"
				 << comp->GetMaterialProperties() << "\t"
				 << endl;
			++number;
		}
	}
}

void ApertureConfiguration::SetLogFile(ostream& os)
{
	log = &os;
}

void ApertureConfiguration::DeleteAllApertures(AcceleratorModel* Model)
{
	std::vector<AcceleratorComponent*> Elements;
	int nElements = Model->ExtractTypedElements(Elements, "*");

	for(const auto &comp : Elements)
	{
		if(comp->GetAperture() != nullptr)
		{
			delete comp->GetAperture();
			comp->SetAperture(nullptr);
		}
	}
}
