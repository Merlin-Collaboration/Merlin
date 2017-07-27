#ifndef DataTableReader_h
#define DataTableReader_h 1

#include <iostream>
#include <fstream>
#include <string>
#include <memory>


#include "DataTable.h"
class DataTable;

class DataTableReader
{
public:
	virtual std::unique_ptr<DataTable> Read() = 0;
private:
};

class DataTableWriter
{
public:
	virtual void Write(DataTable &data_table) = 0;
private:
};

#endif
