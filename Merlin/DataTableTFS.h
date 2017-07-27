#ifndef DataTableTFS_h
#define DataTableTFS_h 1

#include <memory>
#include <iostream>
#include <string>

#include "DataTableReaderWriter.h"
#include "DataTable.h"


class DataTableReaderTFS: public DataTableReader
{
public:
	DataTableReaderTFS(std::istream *in): in(in) {};
	DataTableReaderTFS(std::string filename);

	virtual std::unique_ptr<DataTable> Read() override;

private:
	std::istream *in; // either a passed pointer, or pointer to the opened file
	std::shared_ptr<std::istream> inf; // if we opened the file, this ensures that it is closed
};

class DataTableWriterTFS: public DataTableWriter
{
public:
	DataTableWriterTFS(): width_int(8), width_float(18), prec_float(10) {};
	DataTableWriterTFS(std::ostream *_out):DataTableWriterTFS()
	{
		out = _out;
	};
	DataTableWriterTFS(std::string filename);

	void Write(DataTable& dt);
private:
	int width_int;
	int width_float;
	int prec_float;

	std::ostream *out;

};

#endif
