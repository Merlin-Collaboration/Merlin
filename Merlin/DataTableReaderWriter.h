/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef DataTableReader_h
#define DataTableReader_h 1

#include <iostream>
#include <fstream>
#include <string>
#include <memory>

#include "DataTable.h"
class DataTable;

/** @brief Base class for DataTable readers
 *
 * Derived classes should define a constructor that takes augments needed
 * to set up a source (for example a file name) and a override the Read()
 * method to implement reading.
 */
class DataTableReader
{
public:
	/// Read the source and return a DataTable
	virtual std::unique_ptr<DataTable> Read() = 0;
private:
};

/** @brief Base class for DataTable writers
 *
 * Derived classes should define a constructor that takes augments needed
 * to set up a destination (for example a file name) and a override the
 * Write() method to implement reading.
 */
class DataTableWriter
{
public:
	/// Write the DataTable to destination
	virtual void Write(DataTable &data_table) = 0;
private:
};

#endif
