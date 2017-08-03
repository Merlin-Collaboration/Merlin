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
