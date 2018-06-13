/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef MADKeyMap_h
#define MADKeyMap_h 1

#include "merlin_config.h"
#include <string>
#include <vector>
#include <iostream>
#include <map>

/**
 *      Implementation class for mapping column keys in optics
 *      listing to element types during construction.
 */

class MADKeyMap
{
public:
	virtual ~MADKeyMap()
	{
	}
	typedef std::map<std::string, size_t> key_map;
	struct bad_key {};

	/**
	 * Constructs a key map from the optics listing line
	 * containing the column headings.
	 */
	MADKeyMap(const std::string& hstr);

	/**
	 * Returns the value of the parameter for a specified key.
	 * Throws BAD_KEY if not present.
	 */
	virtual double GetParameter(const std::string& key, bool warn = true);

	/**
	 * Reads in the values for the next row.
	 */
	virtual void ReadRow(std::istream& is);

	bool has_type;
	bool has_apertype;

private:

	std::vector<double> vals;
	std::string type_str;
	key_map kmap;
	size_t apertype_column;
};

#endif
