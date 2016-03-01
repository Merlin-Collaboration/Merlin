#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "AcceleratorModel/Apertures/RectEllipseAperture.h"
#include "AcceleratorModel/Apertures/InterpolatedApertures.h"
#include "AcceleratorModel/StdComponent/Collimator.h"

#include "Collimators/ApertureConfiguration.h"

using namespace std;

ApertureConfiguration::ApertureConfiguration(string input_file) : logFlag(false), allRectEllipse(0)
{
	LoadApertureConfiguration(input_file);
}

ApertureConfiguration::ApertureConfiguration(string input_file, bool are) : logFlag(false), allRectEllipse(are)
{
	LoadApertureConfiguration(input_file);
}


void ApertureConfiguration::LoadApertureConfiguration(string input_file)
{
	ifstream* input = new ifstream(input_file.c_str(), ifstream::in);
	//Do standard checks
	if(input == NULL || !input->good())
	{
		std::cerr << "Failed to open aperture input file: " << input_file << " - Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}
	//string buf,s,ap1,ap2,ap3,ap4;
	string key,name,parent,l,ap1s;
	double s,ap1,ap2,ap3,ap4;

	string aptype;

	//Lets assume we have the correct file layout
	if(input->good())
	{
		char getlinebuf[1024];
		for(size_t n = 47; n > 0; n--)
		{
			input->getline(getlinebuf,1023);
		}
	}

	while(input->good())
	{
		if(!allRectEllipse)
		{
			//				* KEYWORD	NAME	PARENT	S	L	APER_1	APER_2	APER_3	APER_4
			(*input) >> 	key >> 		name >>	parent >> s >> 	l >> 	ap1 >> 	ap2 >> 	ap3 >> 	ap4 >> aptype;
		}
		else
		{
			(*input) >> 	key >> 		name >>	parent >> s >> 	l >> 	ap1 >> 	ap2 >> 	ap3 >> 	ap4;
		}
		if( ( ap1 !=0) || (ap2 !=0) || (ap3 != 0) || (ap4 !=0) )
		{
			ApertureEntry.s = s;
			ApertureEntry.ap1 = ap1;
			ApertureEntry.ap2 = ap2;
			ApertureEntry.ap3 = ap3;
			ApertureEntry.ap4 = ap4;
			if(!allRectEllipse)
			{
				if(aptype == "\"CIRCLE\"")
				{
					//cout << "CIRCLE FIX" << endl;
					ApertureEntry.ap2 = ap1;
					ApertureEntry.ap3 = ap1;
					ApertureEntry.ap4 = ap1;
				}
			}
			ApertureList.push_back(ApertureEntry);
//			cout << ApertureEntry.s << "\t" << ApertureEntry.ap1 << "\t" << ApertureEntry.ap2 << "\t" << ApertureEntry.ap3 << "\t" << ApertureEntry.ap4 << endl;
		}
	}
	delete input;
}

void ApertureConfiguration::ConfigureElementApertures(AcceleratorModel* model)
{
	//Get a list of all elements
	vector<AcceleratorComponent*> elements;
	int nelements = model->ExtractTypedElements(elements,"*");
	std::cout << "Got " << nelements << " elements for aperture configuration" << std::endl;

	/**
	* elements - The container with all AcceleratorComponent entries
	* comp - Iterator to all AcceleratorComponent entries
	*
	* ApertureList - The container with all Aperture entries
	* itr - Iterator to all Aperture entries.
	*
	* ThisElementAperture - The container for the current element Aperture.
	*/

	for(vector<AcceleratorComponent*>::iterator comp = elements.begin(); comp!=elements.end(); comp++)
	{
		//Do not overwrite collimator apertures
		Collimator* collimator = NULL;
		if((*comp)->GetAperture() == NULL)
		{
			double element_length = (*comp)->GetLength();
			double position = (*comp)->GetComponentLatticePosition();

			//We only care about non-zero length elements
			if(element_length != 0)
			{
				//cout <<	(*comp)->GetName() << "\t" << element_length << endl;

				//Loop over all aperture entries
				for(vector<ap>::iterator itr = ApertureList.begin(); itr!=ApertureList.end(); itr++)
				{
					//If the aperture entry is >= the element position
					if(itr->s >= position)
					{
						//gone past where we need to go
						//cout << (*comp)->GetName() << "\t" << itr->s << endl;
						vector<ap> ThisElementAperture;

						//itr-- will give the point before the element + should check for first element if a ring
						//The following does assume symmetry
						if(itr == ApertureList.begin())
						{
							std::cout << "At first element, getting aperture iterpolation from last element" << std::endl;
							//got the initial point
							itr = ApertureList.end();
							itr--;
							(*itr).s = 0;
							//record ap points
							ThisElementAperture.push_back(*itr);
							//go back to where we were
							itr = ApertureList.begin();
							itr++;
						}
						else
						{
							//got the initial point
							itr--;
							//record ap points
							ThisElementAperture.push_back(*itr);
							//go back to where we were
							itr++;
						}


						//Now add in all entries that exist within the length of the element
						while(itr->s <= (position + element_length))
						{
							//Check if we are at the end of the aperture entries - i.e. at the end of the ring.
							if(itr != ApertureList.end())
							{
								//Here we add in a standard element aperture entry and iterate to the next position.
								ThisElementAperture.push_back(*itr);

								//Increment the iterator to the next entry
								itr++;
							}
							else
							{
								//If we are at the last element, then break out of the loop and grab the first entry.
								std::cout << "At last element, getting aperture iterpolation from first element" << std::endl;
								break;
							}
						}

						//grab the last point - for the last element we will need to loop over to the first element again for a ring.
						if(itr == ApertureList.end())
						{
							//cout << "at end: " << (*comp)->GetName() << endl;
							//cout << itr->s << "\t"<< itr->ap2 <<  endl;
							itr = ApertureList.begin();
							(*itr).s += element_length + position;
							std::cout << (*itr).s << std::endl;
							std::cout << (*itr).ap1 << std::endl;
							std::cout << (*itr).ap2 << std::endl;
							std::cout << (*itr).ap3 << std::endl;
							std::cout << (*itr).ap4 << std::endl;

							ThisElementAperture.push_back(*itr);
							//std::cout << "ApertureConfiguration: Failed to obtain interpolation from final element to first" << std::endl;
							//abort();
						}
						else
						{
							ThisElementAperture.push_back(*itr);
						}

						//cout << position << "\t" << (*comp)->GetName()  << "\t" << ThisElementAperture.size() << endl;

						/**
						* Now move on to assigning the correct type of aperture.
						*/

						//aper_# means for all apertypes but racetrack:
						//aper_1 = half width rectangle
						//aper_2 = half heigth rectangle
						//aper_3 = half horizontal axis ellipse (or radius if circle)
						//aper_4 = half vertical axis ellipse

						//Check if all values are constant
						//if so check if we are a circle or rect

						//if not constant check for circle
						//if not circle -> rectellipse
						bool ap1=true,ap2=true,ap3=true,ap4=true,circle=true;

						if(!allRectEllipse)
						{
							circle=true;
						}
						else
						{
							circle = false;
						}

						double ap1p=0,ap2p=0,ap3p=0,ap4p=0;
						//for(vector<ap>::iterator itap = ThisElementAperture.begin(); itap!=ThisElementAperture.end(); itap++)
						for(unsigned int itap = 0; itap < ThisElementAperture.size(); itap++)
						{
							//configure for first pass
							if(itap == 0)
							{
								ap1p = ThisElementAperture[itap].ap1;
								ap2p = ThisElementAperture[itap].ap2;
								ap3p = ThisElementAperture[itap].ap3;
								ap4p = ThisElementAperture[itap].ap4;
							}

							//now check which elements are the same each pass
							if(ThisElementAperture[itap].ap1 != ap1p)
							{
								ap1 = false;
							}
							if(ThisElementAperture[itap].ap2 != ap2p)
							{
								ap2 = false;
							}
							if(ThisElementAperture[itap].ap3 != ap3p)
							{
								ap3 = false;
							}
							if(ThisElementAperture[itap].ap4 != ap4p)
							{
								ap4 = false;
							}
							if(ThisElementAperture[itap].ap1 != ThisElementAperture[itap].ap2)
							{
								circle = false;
							}
							if(ThisElementAperture[itap].ap1 != ThisElementAperture[itap].ap3)
							{
								circle = false;
							}
							if(ThisElementAperture[itap].ap1 != ThisElementAperture[itap].ap4)
							{
								circle = false;
							}

							//always set up for the next pass
							ap1p = ThisElementAperture[itap].ap1;
							ap2p = ThisElementAperture[itap].ap2;
							ap3p = ThisElementAperture[itap].ap3;
							ap4p = ThisElementAperture[itap].ap4;
						}

						//Now to decide what type of aperture to make
						Aperture* aper;
						//cout << (*comp)->GetName() << "\t" << position + element_length << "\t";

						if(circle == true )
						{
							if(ap1 == false || ap2 == false || ap3 == false || ap4 == false)
							{
								//cout << "Interpolated Circle" << endl;

								InterpolatedAperture* apInterpolated = new InterpolatedAperture();

								for(size_t n=0; n < ThisElementAperture.size(); n++ )
								{
									ThisElementAperture[n].s -= position;
									apInterpolated->ApertureEntry.s = ThisElementAperture[n].s;
									apInterpolated->ApertureEntry.ap1 = ThisElementAperture[n].ap1;
									apInterpolated->ApertureEntry.ap2 = ThisElementAperture[n].ap2;
									apInterpolated->ApertureEntry.ap3 = ThisElementAperture[n].ap3;
									apInterpolated->ApertureEntry.ap4 = ThisElementAperture[n].ap4;
									/*
									cout << "ap3: " << ThisElementAperture[n].ap3 << endl;
									cout << "ap3 loaded: " << apInterpolated->ApertureEntry.ap3 << endl;
									cout << "ap1: " << apInterpolated->ApertureEntry.ap1 << endl;
									cout << "ap3: " << apInterpolated->ApertureEntry.ap2 << endl;
									cout << "ap4: " << apInterpolated->ApertureEntry.ap4 << endl << endl;
									*/

									if(ThisElementAperture[n].ap3 == 0 || apInterpolated->ApertureEntry.ap3 == 0)
									{
										for(size_t m=0; m < ThisElementAperture.size(); m++ )
										{
											std::cout << ThisElementAperture[m].s << "\t" << ThisElementAperture[m].ap3 << std::endl;
										}
										abort();
									}
									apInterpolated->ApertureList.push_back(apInterpolated->ApertureEntry);

									if(ThisElementAperture[n].ap4 < 0)
									{
										std::cout << "broken 4" << std::endl;
									}
									//cout << ThisElementAperture[n].s << "\t" << element_length << endl;
								}

								aper = new InterpolatedCircularAperture(apInterpolated->GetApertureList());
								(*comp)->SetAperture(aper);
								ApertureType = "Interpolated Circular";

							}
							else
							{
								//cout << "Circle" << endl;
								aper = new CircularAperture(ap3p);
								(*comp)->SetAperture(aper);
								ApertureType = "Circular";
							}
						}

						else
						{
							if(ap1 == false || ap2 == false || ap3 == false || ap4 == false)
							{
//							cout << "Interpolated RectEllipse" << endl;
								//aper = new InteroplatedRectEllipseAperture(ThisElementAperture);
								//(*comp)->SetAperture(aper);
								InterpolatedAperture* apInterpolated = new InterpolatedAperture();

								for(size_t n=0; n < ThisElementAperture.size(); n++ )
								{
									ThisElementAperture[n].s -= position;
									apInterpolated->ApertureEntry.s = ThisElementAperture[n].s;
									apInterpolated->ApertureEntry.ap1 = ThisElementAperture[n].ap1;
									apInterpolated->ApertureEntry.ap2 = ThisElementAperture[n].ap2;
									apInterpolated->ApertureEntry.ap3 = ThisElementAperture[n].ap3;
									apInterpolated->ApertureEntry.ap4 = ThisElementAperture[n].ap4;
									apInterpolated->ApertureList.push_back(apInterpolated->ApertureEntry);
									if(ThisElementAperture[n].ap4 < 0)
									{
										cout << "broken 4" << endl;
									}
									//cout << ThisElementAperture[n].s << "\t" << element_length << endl;
								}

								aper = new InterpolatedRectEllipseAperture(apInterpolated->GetApertureList());
								(*comp)->SetAperture(aper);
								ApertureType = "Interpolated RectEllipse";
								delete apInterpolated;
							}
							else
							{
								//cout << "RectEllipse" << endl;
								aper = new RectEllipseAperture(ap1p,ap2p,ap3p,ap4p);
								(*comp)->SetAperture(aper);
								ApertureType = "RectEllipse";
							}

						}

						//make aperture class;
						//interpolate for ap1,2,3,4 in turn
						//get ap1
						//get ap1 at s--
						//work out funct
						//will have to create an interpolated aperture that calculates the aperture on the fly
						//give each element the aperture points it requires + 1 either side
						//(*comp)->SetAperture(aper);
						//cout << "Aperture Load end" << endl;
						break;
						delete aper;

					}
				}
			}
		}

		if(logFlag)
		{
			*log << std::setw(25) << std::left << (*comp)->GetName();
			*log << std::setw(14) << std::left << (*comp)->GetType();
			*log << std::setw(10) << std::left << (*comp)->GetLength();
			*log << std::setw(10) << std::left << (*comp)->GetComponentLatticePosition();

			if ((*comp)->GetAperture() != NULL)
			{
				(*comp)->GetAperture()->printout(*log);
			}

			*log << endl;
		}
	}
}


void ApertureConfiguration::SetLogFile (ostream& os)
{
	log=&os;
}

void ApertureConfiguration::EnableLogging(bool flg)
{
	logFlag = flg;
}

