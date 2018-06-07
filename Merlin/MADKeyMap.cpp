/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <sstream>
#include "MerlinIO.h"
#include <cstdlib>

#include "MADKeyMap.h"
using namespace std;

MADKeyMap::MADKeyMap(const std::string& hstr) :
	has_type(false), has_apertype(false)
{
	istringstream is(hstr);
	size_t n = 0;
	string s;

	while(is >> s)
	{
		if(s == "TYPE")
		{
			//This if test is essentially old legacy code and should be removed, but kept to keep compatibility.
			has_type = true;
		}
		else if(s == "APERTYPE")
		{
			//We want to record the column number that contains the aperture type text,
			//for latter conversion from a string to a double.

			//First set that the apertype column exists
			has_apertype = true;

			//Store this value for later usage.
			apertype_column = n;

			//Enter the column number as normal
			kmap[s] = n++;
		}
		else
		{
			kmap[s] = n++;
		}
	}

#ifndef NDEBUG
	cout << n << " column headings identified" << endl;
#endif

	vals = vector<double>(n, 0.0);
}

double MADKeyMap::GetParameter(const std::string& key, bool warn)
{
	key_map::iterator p = kmap.find(key);
	if(p != kmap.end())
	{
		return vals[p->second];
	}

	else
	{
		if(warn)
		{
#ifndef NDEBUG
			MerlinIO::warning() << key << " not in optics listing. Defaulted to zero" << endl;
#endif
		}
		return 0;
	}
}

/*
   Aperture types - from MADX: http://mad.web.cern.ch/mad/Introduction/aperture.html
   CIRCLE		1
   ELLIPSE		2
   RECTANGLE	3
   LHCSCREEN	4
   MARGUERITE	5
   RECTELLIPSE	6
   RACETRACK	7
   NONE		0
 */

void MADKeyMap::ReadRow(std::istream& is)
{
	//Since vals is an array of doubles, we want to convert the aperture text into a
	//number defining the aperture type to be entered into this array
	string buf;
	for(size_t i = 0; i < vals.size(); i++)
	{

		//Are we dealing with the aperture column in the input file?
		if(i == apertype_column && has_apertype == true)
		{
			//~ cout << "MADKeyMap ReadRow() entering has apertype conditional = " << vals[i] << endl;
			is >> buf;
			if(buf == "\"CIRCLE\"")
			{
				//cout << "Have CIRCLE aperture"  << endl;
				vals[i] = 1.0;
			}
			else if(buf == "\"ELLIPSE\"")
			{
				//cout << "Have ELLIPSE aperture"  << endl;
				vals[i] = 2.0;
			}
			else if(buf == "\"RECTANGLE\"")
			{
				//cout << "Have RECTANGLE aperture"  << endl;
				vals[i] = 3.0;
			}
			else if(buf == "\"LHCSCREEN\"")
			{
				//cout << "Have LHCSCREEN aperture"  << endl;
				vals[i] = 4.0;
			}
			else if(buf == "\"MARGUERITE\"")
			{
				//cout << "Have MARGUERITE aperture"  << endl;
				vals[i] = 5.0;
			}
			else if(buf == "\"RECTELLIPSE\"")
			{
				//cout << "Have RECTELLIPSE aperture"  << endl;
				vals[i] = 6.0;
			}
			else if(buf == "\"RACETRACK\"")
			{
				//cout << "Have RACETRACK aperture"  << endl;
				vals[i] = 7.0;
			}
			else if(buf == "\"NONE\"")
			{
				//cout << "Have no aperture!"  << endl;
				vals[i] = 0.0;
			}

		}
		//else read as normal
		else
		{
			is >> vals[i];
		}
	}
}
