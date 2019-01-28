/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "DataTable.h"
#include <iostream>

void DataTable::AddColumn(std::string col_name, char type)
{
	if(lookup.count(col_name))
	{
		throw std::invalid_argument("Column '" + col_name + "' already exists");
	}

	size_t position;
	switch(type)
	{
	case 'd':
	{
		position = data_d.size();
		data_d.emplace_back(length);
		break;
	}
	case 'i':
	{
		position = data_i.size();
		data_i.emplace_back(length);
		break;
	}
	case 's':
	{
		position = data_s.size();
		data_s.emplace_back(length);
		break;
	}
	default:
		throw WrongTypeException(std::string("Unknown type: '") + type + "'");
	}

	col_names.push_back(col_name);
	lookup[col_name].type = type;
	lookup[col_name].pos = position;
}

size_t DataTable::AddRow()
{
	for(auto &a : data_d)
	{
		a.emplace_back();
	}
	for(auto &a : data_i)
	{
		a.emplace_back();
	}
	for(auto &a : data_s)
	{
		a.emplace_back();
	}
	length++;
	return length - 1;
}

double DataTable::Get_d(const std::string col_name, size_t i) const
{
	location l = lookup.at(col_name);
	if(l.type != 'd')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a double");
	}
	return data_d.at(l.pos).at(i);
}

int DataTable::Get_i(const std::string col_name, size_t i) const
{
	location l = lookup.at(col_name);
	if(l.type != 'i')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a int");
	}
	return data_i.at(l.pos).at(i);
}

std::string DataTable::Get_s(const std::string col_name, size_t i) const
{
	location l = lookup.at(col_name);
	if(l.type != 's')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a string");
	}
	return data_s.at(l.pos).at(i);
}

void DataTable::Set(const std::string col_name, size_t i, double x)
{
	location l = lookup.at(col_name);
	if(l.type != 'd')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a double");
	}
	data_d.at(l.pos).at(i) = x;
}

void DataTable::Set(const std::string col_name, size_t i, int x)
{
	location l = lookup.at(col_name);
	if(l.type != 'i')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a int");
	}
	data_i.at(l.pos).at(i) = x;
}

void DataTable::Set(const std::string col_name, size_t i, std::string x)
{
	location l = lookup.at(col_name);
	if(l.type != 's')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a string");
	}
	data_s.at(l.pos).at(i) = x;
}

void DataTable::SetWithStr(const std::string col_name, size_t i, std::string x)
{
	location l = lookup.at(col_name);
	switch(l.type)
	{
	case 'd':
		Set_d(col_name, i, stod(x));
		break;
	case 'i':
		Set_i(col_name, i, stoi(x));
		break;
	case 's':
		Set_s(col_name, i, x);
		break;
	default:
		throw WrongTypeException("Bad type in lookup");
	}
}

std::string DataTable::GetAsStr(const std::string col_name, size_t i) const
{
	location l = lookup.at(col_name);
	switch(l.type)
	{
	case 'd':
		return std::to_string(Get_d(col_name, i));
	case 'i':
		return std::to_string(Get_i(col_name, i));
	case 's':
		return Get_s(col_name, i);
	default:
		throw WrongTypeException("Bad type in lookup");
	}
}

void DataTable::HeaderAddColumn(std::string col_name, char type)
{
	if(hlookup.count(col_name))
	{
		throw std::invalid_argument("Header '" + col_name + "' already exists");
	}

	size_t position;
	switch(type)
	{
	case 'd':
	{
		position = hdata_d.size();
		hdata_d.emplace_back();
		break;
	}
	case 'i':
	{
		position = hdata_i.size();
		hdata_i.emplace_back();
		break;
	}
	case 's':
	{
		position = hdata_s.size();
		hdata_s.emplace_back();
		break;
	}
	default:
		throw WrongTypeException(std::string("Unknown type: '") + type + "'");
	}

	hcol_names.push_back(col_name);
	hlookup[col_name].type = type;
	hlookup[col_name].pos = position;
}

double DataTable::HeaderGet_d(const std::string col_name) const
{
	location l = hlookup.at(col_name);
	if(l.type != 'd')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a double");
	}
	return hdata_d.at(l.pos);
}

int DataTable::HeaderGet_i(const std::string col_name) const
{
	location l = hlookup.at(col_name);
	if(l.type != 'i')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a int");
	}
	return hdata_i.at(l.pos);
}

std::string DataTable::HeaderGet_s(const std::string col_name) const
{
	location l = hlookup.at(col_name);
	if(l.type != 's')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a string");
	}
	return hdata_s.at(l.pos);
}

void DataTable::HeaderSet(const std::string col_name, double x)
{
	location l = hlookup.at(col_name);
	if(l.type != 'd')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a double");
	}
	hdata_d.at(l.pos) = x;
}

void DataTable::HeaderSet(const std::string col_name, int x)
{
	location l = hlookup.at(col_name);
	if(l.type != 'i')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a int");
	}
	hdata_i.at(l.pos) = x;
}

void DataTable::HeaderSet(const std::string col_name, std::string x)
{
	location l = hlookup.at(col_name);
	if(l.type != 's')
	{
		throw WrongTypeException("Column '" + col_name + "' is not a string");
	}
	hdata_s.at(l.pos) = x;
}

void DataTable::HeaderSetWithStr(const std::string col_name, std::string x)
{
	location l = hlookup.at(col_name);
	switch(l.type)
	{
	case 'd':
		HeaderSet_d(col_name, stod(x));
		break;
	case 'i':
		HeaderSet_i(col_name, stoi(x));
		break;
	case 's':
		HeaderSet_s(col_name, x);
		break;
	default:
		throw WrongTypeException("Bad type in lookup");
	}
}

std::string DataTable::HeaderGetAsStr(const std::string col_name) const
{
	location l = hlookup.at(col_name);
	switch(l.type)
	{
	case 'd':
		return std::to_string(HeaderGet_d(col_name));
	case 'i':
		return std::to_string(HeaderGet_i(col_name));
	case 's':
		return HeaderGet_s(col_name);
	default:
		throw WrongTypeException("Bad type in lookup");
	}
}

bool DataTable::HasCol(std::string col_name) const
{
	return lookup.count(col_name) == 1;
}

bool DataTable::HasCol(std::string col_name, char type) const
{
	return lookup.count(col_name) == 1 && lookup.at(col_name).type == type;
}

bool DataTable::HeaderHasKey(std::string col_name) const
{
	return hlookup.count(col_name) == 1;
}

bool DataTable::HeaderHasKey(std::string col_name, char type) const
{
	return hlookup.count(col_name) == 1 && hlookup.at(col_name).type == type;
}

void DataTable::OutputAscii(std::ostream &os) const
{
	for(auto &col_name : hcol_names)
	{
		os << col_name << "\t" << HeaderGetAsStr(col_name) << std::endl;
	}

	for(auto &col_name : col_names)
	{
		os << col_name << "\t";
	}
	os << std::endl;
	for(size_t i = 0; i < length; i++)
	{
		for(auto &col_name : col_names)
		{
			os << GetAsStr(col_name, i) << "\t";
		}
		os << std::endl;
	}
}

DataTableRowIterator DataTable::begin() const
{
	return DataTableRowIterator(this, 0);
}

DataTableRowIterator DataTable::end() const
{
	return DataTableRowIterator(this, length);
}

DataTableRow DataTableRowIterator::operator *()
{
	return DataTableRow(dt, pos);
}

DataTableRowPtr DataTableRowIterator::operator ->()
{
	return DataTableRowPtr(dt, pos);
}

DataTableRowIterator& DataTableRowIterator::operator++()
{
	++pos;
	return *this;
}

DataTableRowIterator DataTableRowIterator::operator++(int)
{
	auto ret = *this;
	++pos;
	return ret;
}

DataTableRowIterator& DataTableRowIterator::operator--()
{
	--pos;
	return *this;
}

DataTableRowIterator DataTableRowIterator::operator--(int)
{
	auto ret = *this;
	--pos;
	return ret;
}

bool DataTableRowIterator::operator==(const DataTableRowIterator &other) const
{
	return (dt == other.dt) && (pos == other.pos);
}

bool DataTableRowIterator::operator!=(const DataTableRowIterator &other) const
{
	return !(*this == other);
}
