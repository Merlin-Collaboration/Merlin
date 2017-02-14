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

ApertureConfiguration::ApertureConfiguration(std::string InputFileName) : logFlag(false), DefaultAperture(nullptr), DefaultApertureFlag(false)

{
	LoadApertureConfiguration(InputFileName);
}
/*
ApertureConfiguration::ApertureConfiguration(std::string InputFileName, bool are) : logFlag(false), allRectEllipse(are)
{
	LoadApertureConfiguration(InputFileName);
}
*/

void ApertureConfiguration::LoadApertureConfiguration(std::string InputFileName)
{
	std::ifstream* input = new std::ifstream(InputFileName.c_str(), std::ifstream::in);

	//Do standard checks
	if(input == nullptr || !input->good())
	{
		std::cerr << "Failed to open aperture input file: " << InputFileName << " - Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string key,name,parent,aptype;
	double s,l,ap1,ap2,ap3,ap4;

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
		//				* KEYWORD	NAME	PARENT	S	L	APER_1	APER_2	APER_3	APER_4  APTYPE
		(*input) >> 	key >> 		name >>	parent >> s >> 	l >> 	ap1 >> 	ap2 >> 	ap3 >> 	ap4 >> aptype;

		if( ( ap1 !=0) || (ap2 !=0) || (ap3 != 0) || (ap4 !=0) )
		{
			//Assume refer = EXIT (if not we are in trouble)
			ApertureEntry.s = s-l;
			ApertureEntry.ap1 = ap1;
			ApertureEntry.ap2 = ap2;
			ApertureEntry.ap3 = ap3;
			ApertureEntry.ap4 = ap4;
			ApertureEntry.ap4 = ap4;

			if(aptype == "\"RECTANGLE\"")
			{
				ApertureEntry.ApType = RECTANGLE;
			}
			else if(aptype == "\"CIRCLE\"")
			{
				ApertureEntry.ApType = CIRCLE;
			}
			else if(aptype == "\"ELLIPSE\"")
			{
				ApertureEntry.ApType = ELLIPSE;
			}
			else if(aptype == "\"RECTELLIPSE\"")
			{
				ApertureEntry.ApType = RECTELLIPSE;
			}
			else if(aptype == "\"LHCSCREEN\"")
			{
				ApertureEntry.ApType = LHCSCREEN;
			}
			else if(aptype == "\"OCTAGON\"")
			{
				ApertureEntry.ApType = OCTAGON;
			}
			else if(aptype == "\"NONE\"")
			{
				ApertureEntry.ApType = NONE;
			}
			else
			{
				//Unknown type
				std::cerr << "Unknown Aperture type found in Aperture input file: " << aptype << " - Exiting!" << std::endl;
				exit(EXIT_FAILURE);
			}

			//Front of element
			ApertureList.push_back(ApertureEntry);

			//Exit of element
			ApertureEntry.s = s;
			ApertureList.push_back(ApertureEntry);
		}
	}
	if(input)
	{
		input->close();
	}
	delete input;
}

void ApertureConfiguration::OutputApertureList(std::ostream& os)
{
	for(size_t n=0; n < ApertureList.size(); n++)
	{
		os << ApertureList[n].s << "\t" << ApertureList[n].ap1 << "\t" << ApertureList[n].ap2 << "\t" << ApertureList[n].ap3 << "\t" << ApertureList[n].ap4 << "\t";
		if(ApertureList[n].ApType == RECTELLIPSE)
		{
			os << "RECTELLIPSE";
		}
		else if(ApertureList[n].ApType == CIRCLE)
		{
			os << "CIRCLE";
		}
		else if(ApertureList[n].ApType == RECTANGLE)
		{
			os << "RECTANGLE";
		}
		else if(ApertureList[n].ApType == LHCSCREEN)
		{
			os << "LHCSCREEN";
		}
		else if(ApertureList[n].ApType == ELLIPSE)
		{
			os << "ELLIPSE";
		}
		else if(ApertureList[n].ApType == OCTAGON)
		{
			os << "OCTAGON";
		}
		else
		{
			os << "BUG";
		}
		os << std::endl;
	}
}

void ApertureConfiguration::ConfigureElementApertures(AcceleratorModel* Model)
{
	//Get a list of all elements
	std::vector<AcceleratorComponent*> Elements;
	int nElements = Model->ExtractTypedElements(Elements,"*");
	std::cout << "Got " << nElements << " elements for aperture configuration" << std::endl;
	std::cout << "Got " << ApertureList.size() << " Aperture entries" << std::endl;

	Aperture* aper;
	/**
	* Elements - The container with all AcceleratorComponent entries
	* comp - Iterator to all AcceleratorComponent entries
	*
	* ApertureList - The container with all Aperture entries
	* itr - Iterator to all Aperture entries.
	*
	* ThisElementAperture - The container for the current element Aperture.
	*/

	for(std::vector<AcceleratorComponent*>::iterator comp = Elements.begin(); comp!=Elements.end(); comp++)
	{
		if((*comp)->GetAperture() == nullptr)
		{
			double ElementLength = (*comp)->GetLength();
			double Position = (*comp)->GetComponentLatticePosition();

			/**
			* Give magnets and other fixed elements a set aperture.
			* Give drifts interpolated apertures.
			*/

			//Get the Type of Element this is.
			std::string ElementType = (*comp)->GetType();

			//We only care about non-zero length elements
			if(ElementLength != 0)
			{
				//For all elements that are not drifts, set a fixed (non-interpolated) aperture type.
				if(ElementType != "Drift")
				{
					//std::cout << (*comp)->GetQualifiedName() << " run" << std::endl;
					//Loop over all aperture entries
					for(std::vector<ap>::iterator itr = ApertureList.begin(); itr!=ApertureList.end(); itr++)
					{
						//If the aperture entry is > the element start position and less than the exit
						//if(( itr->s >= Position && itr->s <= (Position + ElementLength) ) || fequal(itr->s, (Position + ElementLength),1e-6) )
						//if(( itr->s >= Position)  || fequal(itr->s, (Position + ElementLength),1e-6) )
						if( itr->s >= Position )
						{
							//	std::cout << (*comp)->GetQualifiedName() << " hit" << std::endl;
							//Add the appropriate aperture type
							if(itr->ApType == RECTELLIPSE || itr->ApType == LHCSCREEN)
							{
								aper = new RectEllipseAperture(itr->ap1, itr->ap2, itr->ap3, itr->ap4);
							}
							else if(itr->ApType == CIRCLE)
							{
								aper = new CircularAperture(itr->ap1);
							}
							else if(itr->ApType == ELLIPSE)
							{
								aper = new EllipticalAperture(itr->ap1, itr->ap2);
							}
							else if(itr->ApType == RECTANGLE)
							{
								aper = new RectangularAperture(itr->ap1, itr->ap2);
							}
							else if(itr->ApType == OCTAGON)
							{
								aper = new OctagonalAperture(itr->ap1, itr->ap2, itr->ap3, itr->ap4);
							}
							else
							{
								std::cerr << (*comp)->GetQualifiedName() << " aperture Class bug" << std::endl;
								exit(EXIT_FAILURE);
							}
							(*comp)->SetAperture(aper);
							itr=ApertureList.end();
							break;
						}
						if(itr==ApertureList.end())
						{
							std::cerr << (*comp)->GetQualifiedName() << " hit end" << std::endl;
						}
					}
				}//End of fixed elements
				else //Deal with drifts
				{
					//std::cout << (*comp)->GetQualifiedName() << std::endl;
					//Loop over all aperture entries
					for(std::vector<ap>::iterator itr = ApertureList.begin(); itr!=ApertureList.end(); itr++)
					{
						//If the aperture entry is >= the element position
						//or if we have reach the end of the aperture entries i.e. the left overs at the end of the ring past the last marker
						if(itr->s >= Position || itr==(ApertureList.end()-1))
						{
							/**
							* 3 possible cases
							* 1: First entry in the element (or accelerator)
							* 2: Entries within an element
							* 3: Final entry within an element (or accelerator)
							*/
							std::vector<ap> ThisElementAperture;

							//Deal with the first entry for this element
							if(itr == ApertureList.begin())
							{
								std::cout << "At first element " << (*comp)->GetQualifiedName() << " getting aperture iterpolation from last element" << std::endl;

								//got the initial point
								itr = ApertureList.end()-1;

								ap tempAp;
								tempAp.s = 0;
								tempAp.ap1 = itr->ap1;
								tempAp.ap2 = itr->ap2;
								tempAp.ap3 = itr->ap3;
								tempAp.ap4 = itr->ap4;
								tempAp.ApType = itr->ApType;

								//record ap points
								ThisElementAperture.push_back(tempAp);

								//go back to where we were
								itr = ApertureList.begin();
							}
							else
							{
								//get the previous point before this element
								itr--;

								//record ap points
								ThisElementAperture.push_back(*itr);

								//go back to where we were
								itr++;
							}

							//Now add in all entries that exist within the length of the element
							while(itr->s <= (Position + ElementLength) )
								//while(itr->s <= (Position + ElementLength))
							{
								ThisElementAperture.push_back(*itr);

								//Increment the iterator to the next entry
								itr++;

								if(itr == ApertureList.end())
								{
									break;
								}
							}

							//Deal with the very last element entry
							if(itr == ApertureList.end())
							{
								std::cout << "At last element " << (*comp)->GetQualifiedName() << " getting aperture iterpolation from first element" << std::endl;
								itr = ApertureList.begin();

								ap tempAp;
								tempAp.s = ElementLength + Position;
								tempAp.ap1 = itr->ap1;
								tempAp.ap2 = itr->ap2;
								tempAp.ap3 = itr->ap3;
								tempAp.ap4 = itr->ap4;
								tempAp.ApType = itr->ApType;

								ThisElementAperture.push_back(tempAp);
							}
							else
							{
								ThisElementAperture.push_back(*itr);
								itr = ApertureList.end();
							}

							/**
							* Now move on to assigning the correct type of aperture.
							*/

							bool ZeroEntry = false;
							size_t NegativeCount = 0;
							//First do a little bit of cleaning
							//If we have an entry at 0 (or very close to), and also an entry at negative values, we can discard the negative entry
							for(size_t itAp = 0; itAp < ThisElementAperture.size(); itAp++)
							{
								if( fequal(ThisElementAperture[itAp].s - Position, 0.0, 1e-7) )
								{
									ZeroEntry = true;
									ThisElementAperture[itAp].s = Position;
								}

								if( ThisElementAperture[itAp].s - Position < 0 )
								{
									NegativeCount++;
								}
							}

							if(NegativeCount!=0 && ZeroEntry)
							{
								//Delete the first entry (negative)
								ThisElementAperture.erase(ThisElementAperture.begin());
								NegativeCount--;
							}

							while(NegativeCount > 1 )
							{
								//Delete the first entry (negative)
								ThisElementAperture.erase(ThisElementAperture.begin());
								NegativeCount--;
							}

							if( fequal(ThisElementAperture[0].s - Position, 0.0, 5e-7) )
							{
								ThisElementAperture[0].s = Position;
							}

							/**
							* Check if all values are constant
							*/
							bool ap1=true,ap2=true,ap3=true,ap4=true,circle=true;
							double ap1p=0,ap2p=0,ap3p=0,ap4p=0;
							bool ApTypeChange = false;
							ApertureClass_t ApTypeToAdd;

							for(size_t itap = 0; itap < ThisElementAperture.size(); itap++)
							{
								//configure for first pass
								if(itap == 0)
								{
									ap1p = ThisElementAperture[itap].ap1;
									ap2p = ThisElementAperture[itap].ap2;
									ap3p = ThisElementAperture[itap].ap3;
									ap4p = ThisElementAperture[itap].ap4;
									ApTypeToAdd = ThisElementAperture[itap].ApType;
								}
								else
								{
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
									if(ThisElementAperture[itap].ApType != ApTypeToAdd)
									{
										ApTypeChange = true;
									}
								}
							}

							bool Interpolated = false;
							if(ap1 == false || ap2 == false || ap3 == false || ap4 == false)
							{
								Interpolated = true;
							}

							if(Interpolated)
							{
								InterpolatedAperture* apInterpolated = new InterpolatedAperture();

								for(size_t n=0; n < ThisElementAperture.size(); n++ )
								{
									//ThisElementAperture[n].s -= Position;
									apInterpolated->ApertureEntry.s = ThisElementAperture[n].s - Position;
									apInterpolated->ApertureEntry.ap1 = ThisElementAperture[n].ap1;
									apInterpolated->ApertureEntry.ap2 = ThisElementAperture[n].ap2;
									apInterpolated->ApertureEntry.ap3 = ThisElementAperture[n].ap3;
									apInterpolated->ApertureEntry.ap4 = ThisElementAperture[n].ap4;

									apInterpolated->ApertureList.push_back(apInterpolated->ApertureEntry);
								}

								//Check we have a constant aperture type
								if(ApTypeChange == false)
								{
									if(ApTypeToAdd == CIRCLE)
									{
										aper = new InterpolatedCircularAperture(apInterpolated->GetApertureList());
									}
									else if(ApTypeToAdd == RECTANGLE)
									{
										std::cerr << "TODO: Add a InterpolatedRectangularAperture Class" << std::endl;
										exit(EXIT_FAILURE);
//										aper = new InterpolatedRectangularAperture(apInterpolated->GetApertureList());
									}
									else if(ApTypeToAdd == ELLIPSE)
									{
										aper = new InterpolatedEllipticalAperture(apInterpolated->GetApertureList());
									}
									else if(ApTypeToAdd == RECTELLIPSE || ApTypeToAdd == LHCSCREEN)
									{
										aper = new InterpolatedRectEllipseAperture(apInterpolated->GetApertureList());
									}
									else if(ApTypeToAdd == OCTAGON)
									{
										aper = new InterpolatedOctagonalAperture(apInterpolated->GetApertureList());
									}
									else
									{
										std::cerr << "Drift: constant type aperture Class bug: " << (*comp)->GetQualifiedName() << " at " << (*comp)->GetComponentLatticePosition() << "m" << std::endl;
										std::cerr << "Trying to make an INTERPOLATED APERTURE class of a type that does not exist currently!" << std::endl;
										exit(EXIT_FAILURE);
									}
									(*comp)->SetAperture(aper);
									itr = ApertureList.end();
									break;
								}
								else //We have a change in aperture type. Assume rectellipse for now
								{
									//This should work for changes between circles/ellipses/rectellipse
									//Just set the missing coordinates to the same size as the known parameters for rectellipse
									InterpolatedAperture* apInterpolated = new InterpolatedAperture();

									bool octagon = false;
									bool rectellipse = false;

									for(size_t n=0; n < ThisElementAperture.size(); n++ )
									{
										ThisElementAperture[n].s -= Position;

										apInterpolated->ApertureEntry.s = ThisElementAperture[n].s;

										if(ThisElementAperture[n].ApType == RECTELLIPSE)
										{
											apInterpolated->ApertureEntry.ap1 = ThisElementAperture[n].ap1;
											apInterpolated->ApertureEntry.ap2 = ThisElementAperture[n].ap2;
											apInterpolated->ApertureEntry.ap3 = ThisElementAperture[n].ap3;
											apInterpolated->ApertureEntry.ap4 = ThisElementAperture[n].ap4;
											rectellipse = true;
										}
										else if(ThisElementAperture[n].ApType == CIRCLE)
										{
											apInterpolated->ApertureEntry.ap1 = ThisElementAperture[n].ap1;
											apInterpolated->ApertureEntry.ap2 = ThisElementAperture[n].ap1;
											apInterpolated->ApertureEntry.ap3 = ThisElementAperture[n].ap1;
											apInterpolated->ApertureEntry.ap4 = ThisElementAperture[n].ap1;
											rectellipse = true;
										}
										else if(ThisElementAperture[n].ApType == ELLIPSE)
										{
											apInterpolated->ApertureEntry.ap1 = ThisElementAperture[n].ap1;
											apInterpolated->ApertureEntry.ap2 = ThisElementAperture[n].ap2;
											apInterpolated->ApertureEntry.ap3 = ThisElementAperture[n].ap1;
											apInterpolated->ApertureEntry.ap4 = ThisElementAperture[n].ap2;
											rectellipse = true;
										}
										else if(ThisElementAperture[n].ApType == RECTANGLE)
										{
											apInterpolated->ApertureEntry.ap1 = ThisElementAperture[n].ap1;
											apInterpolated->ApertureEntry.ap2 = ThisElementAperture[n].ap2;
											apInterpolated->ApertureEntry.ap3 = ThisElementAperture[n].ap1;
											apInterpolated->ApertureEntry.ap4 = ThisElementAperture[n].ap2;
											rectellipse = true;
										}
										else if(ThisElementAperture[n].ApType == OCTAGON)
										{
											apInterpolated->ApertureEntry.ap1 = ThisElementAperture[n].ap1;
											apInterpolated->ApertureEntry.ap2 = ThisElementAperture[n].ap2;
											apInterpolated->ApertureEntry.ap3 = ThisElementAperture[n].ap3;
											apInterpolated->ApertureEntry.ap4 = ThisElementAperture[n].ap4;
											octagon = true;
										}
										else
										{
											std::cerr << "Drift: changing type interpolated aperture Class bug: " << (*comp)->GetQualifiedName() << " at " << (*comp)->GetComponentLatticePosition() << "m" << std::endl;
											std::cerr << "Trying to make an INTERPOLATED APERTURE class of a type that does not exist currently!: " << ThisElementAperture[n].ApType << std::endl;
											exit(EXIT_FAILURE);
										}

										if(rectellipse && octagon)
										{
											if(DefaultApertureFlag && DefaultAperture)
											{
												aper = DefaultAperture;
											}
											else
											{
												std::cerr << "Drift: changing type interpolated aperture Class bug: " << (*comp)->GetQualifiedName() << " at " << (*comp)->GetComponentLatticePosition() << "m" << std::endl;
												std::cerr << "Trying to connect octagon apertures with types that are not compatible and no default aperture class is set." << std::endl;
												exit(EXIT_FAILURE);
											}
										}

										apInterpolated->ApertureList.push_back(apInterpolated->ApertureEntry);
									}

									aper = new InterpolatedRectEllipseAperture(apInterpolated->GetApertureList());
									(*comp)->SetAperture(aper);
									itr = ApertureList.end();
									break;
								}
							}
							else //constant aperture drift, operate as for magnets
							{
								if(ApTypeToAdd == RECTELLIPSE || ApTypeToAdd == LHCSCREEN)
								{
									aper = new RectEllipseAperture(ap1p, ap2p, ap3p, ap4p);
								}
								else if(ApTypeToAdd == CIRCLE)
								{
									aper = new CircularAperture(ap1p);
								}
								else if(ApTypeToAdd == ELLIPSE)
								{
									aper = new EllipticalAperture(ap1p, ap2p);
								}
								else if(ApTypeToAdd == RECTANGLE)
								{
									aper = new RectangularAperture(ap1p, ap2p);
								}
								else if(ApTypeToAdd == OCTAGON)
								{
									aper = new OctagonalAperture(ap1p, ap2p, ap3p, ap4p);
								}
								else
								{
									std::cerr << "Trying to create an unknown aperture type! ApTypeToAdd = "<< ApTypeToAdd << std::endl;
									exit(EXIT_FAILURE);
								}
								(*comp)->SetAperture(aper);
								itr = ApertureList.end();
								break;
							}
						}
					}

				}//End of drifts
			}
		}

		if(logFlag)
		{
			*log << std::setw(25) << std::left << (*comp)->GetName();
			*log << std::setw(14) << std::left << (*comp)->GetType();
			*log << std::setw(10) << std::left << (*comp)->GetLength();
			*log << std::setw(10) << std::left << (*comp)->GetComponentLatticePosition();

			if ((*comp)->GetAperture() != nullptr)
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

void ApertureConfiguration::DeleteAllApertures(AcceleratorModel* Model)
{
	std::vector<AcceleratorComponent*> Elements;
	int nElements = Model->ExtractTypedElements(Elements,"*");

	for(std::vector<AcceleratorComponent*>::iterator comp = Elements.begin(); comp!=Elements.end(); comp++)
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
