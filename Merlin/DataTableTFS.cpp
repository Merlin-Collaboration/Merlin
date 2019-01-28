/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "DataTableTFS.h"

DataTableReaderTFS::DataTableReaderTFS(std::string filename)
{
	inf = std::make_shared<std::ifstream>(filename);
	if(!inf->good())
	{
		std::cerr << "Could not open file " << filename << std::endl;
		exit(1);
	}
	in = inf.get();
}

static char type_conv(std::string s)
{
	if(s == "%le")
	{
		return 'd';
	}
	else if(s == "%d")
	{
		return 'i';
	}
	else if(s == "%s")
	{
		return 's';
	}
	else if(s.size() > 2 && s[s.size() - 1] == 's')
	{
		return 's';
	}
	throw BadFormatException("Unknown data type'" + s + "'");
}

static const std::string whitespace = " \t\f\v\n\r";

static std::vector<std::string> split_line(std::string line)
{
	// Read line into vector of strings
	// split on whitespace unless quoted
	std::vector<std::string> sl;
	std::string current;
	current.reserve(64);
	bool double_quote_on = false;

	for(auto c : line)
	{
		if(c == '"')
		{
			double_quote_on = !double_quote_on;
		}
		else if(!double_quote_on && whitespace.find(c) != std::string::npos)
		{
			if(current.length())
			{
				sl.push_back(current);
				current = "";
			}
		}
		else
		{
			current += c;
		}
	}
	if(current.length())
	{
		sl.push_back(current);
	}

	return sl;
}

std::unique_ptr<DataTable> DataTableReaderTFS::Read()
{
	std::unique_ptr<DataTable> dt(new DataTable);
	std::vector<std::string> col_names;
	std::vector<char> col_types;
	std::string line;
	std::vector<std::string> words;

	// Read header
	// Lines start with "@". Ends when line starts with "*"
	while(getline(*in, line))
	{
		words = split_line(line);
		if(words.size() == 0)
		{
			continue;
		}
		if(words[0] != "@")
		{
			break;
		}
		if(words.size() != 4)
		{
			throw BadFormatException("Expected '@ name %type value'");
		}

		dt->HeaderAddColumn(words[1], type_conv(words[2]));
		dt->HeaderSetWithStr(words[1], words[3]);
	}

	//Read column names
	if(words[0] != "*")
	{
		throw BadFormatException("Expected line starting with '*' or '@'");
	}
	for(auto w = words.begin() + 1; w < words.end(); w++)
	{
		col_names.push_back(*w);
	}

	//Read column types
	getline(*in, line);
	words = split_line(line);
	if(words[0] != "$")
	{
		throw BadFormatException("Expected line starting with '$'");
	}
	for(auto w = words.begin() + 1; w != words.end(); w++)
	{
		col_types.push_back(type_conv(*w));
	}

	if(col_names.size() != col_types.size())
	{
		throw BadFormatException("Mismatched length of column names and types");
	}

	{
		auto name = col_names.begin();
		auto type = col_types.begin();
		while(name != col_names.end())
		{
			dt->AddColumn(*name, *type);
			++name;
			++type;
		}
	}

	// Read body
	while(getline(*in, line))
	{
		words = split_line(line);
		if(words.size() == 0)
		{
			continue;
		}

		if(col_names.size() != col_types.size())
		{
			throw BadFormatException("Row does not contain correct number of values");
		}
		size_t rown = dt->AddRow();
		{
			auto name = col_names.begin();
			auto word = words.begin();
			while(name != col_names.end())
			{
				dt->SetWithStr(*name, rown, *word);
				++name;
				++word;
			}
		}
	}

	return dt;
}

DataTableWriterTFS::DataTableWriterTFS(std::string filename) :
	DataTableWriterTFS()
{
	outf = std::make_shared<std::ofstream>(filename);
	if(!outf->good())
	{
		std::cerr << "Could not open file " << filename << std::endl;
		exit(1);
	}
	out = outf.get();
}

void DataTableWriterTFS::Write(DataTable & dt)
{
	std::ios init(NULL);
	init.copyfmt(*out);
	*out << std::setprecision(prec_float);

	for(auto &header_name : dt.HeaderNames())
	{
		*out << "@ " << std::setw(16) << std::left << header_name << " ";
		char head_type = dt.GetHeaderType(header_name);
		switch(head_type)
		{
		case 'd':
			*out << "%le ";
			*out << std::setw(width_float) << dt.HeaderGet_d(header_name) << std::endl;
			break;
		case 'i':
			*out << "%d ";
			*out << std::setw(width_int) << dt.HeaderGet_i(header_name) << std::endl;
			break;
		case 's':
			*out << "%" << dt.HeaderGet_s(header_name).size() << "s ";
			*out << '"' << dt.HeaderGet_s(header_name) << '"' << std::endl;
			break;
		}
	}

	*out << "* " << std::left;
	for(auto &col_name : dt.ColumnNames())
	{
		*out << col_name << " ";
	}
	*out << std::endl;

	*out << "$ ";
	for(auto &col_name : dt.ColumnNames())
	{
		char col_type = dt.GetColumnType(col_name);
		switch(col_type)
		{
		case 'd':
			*out << "%le ";
			break;
		case 'i':
			*out << "%d ";
			break;
		case 's':
			*out << "%s ";
			break;

		}
	}
	*out << std::endl;

	auto col_names = dt.ColumnNames();
	for(auto drow : dt)
	{
		for(auto &col_name : col_names)
		{
			char col_type = dt.GetColumnType(col_name);
			*out << " ";
			switch(col_type)
			{
			case 'd':
				*out << std::setw(width_float) << drow.Get_d(col_name);
				break;
			case 'i':
				*out << std::setw(width_int) << drow.Get_i(col_name);
				break;
			case 's':
				*out << '"' << drow.Get_s(col_name) << '"';
				break;
			}
		}
		*out << std::endl;
	}

	out->copyfmt(init);
}
