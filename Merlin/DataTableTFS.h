/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef DataTableTFS_h
#define DataTableTFS_h 1

#include <memory>
#include <iostream>
#include <string>

#include "DataTableReaderWriter.h"
#include "DataTable.h"

/** @brief Read a DataTable from a TFS file
 *
 * For example to read file generated with MadX
 */
class DataTableReaderTFS: public DataTableReader
{
public:
	/// Read from an istream, e.g. an already opened file
	DataTableReaderTFS(std::istream *in) :
		in(in)
	{
	}
	/// Open a file to read
	DataTableReaderTFS(std::string filename);

	/// Read the file, returning a new DataTable
	virtual std::unique_ptr<DataTable> Read() override;

private:
	std::istream *in; // either a passed pointer, or pointer to the opened file
	std::shared_ptr<std::istream> inf; // if we opened the file, this ensures that it is closed
};

/** @brief Write a DataTable to a TFS file
 */
class DataTableWriterTFS: public DataTableWriter
{
public:
	/// Write to an ostream, e.g. an already opened file
	DataTableWriterTFS(std::ostream *_out) :
		DataTableWriterTFS()
	{
		out = _out;
	}
	/// Open a file to write to
	DataTableWriterTFS(std::string filename);

	/// Write the DataTable to the file or stream
	void Write(DataTable& dt);
private:
	DataTableWriterTFS() :
		width_int(8), width_float(18), prec_float(10)
	{
	}
	int width_int;
	int width_float;
	int prec_float;

	std::ostream *out;
	std::shared_ptr<std::ostream> outf;
};

#endif
