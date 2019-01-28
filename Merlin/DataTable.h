/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef DataTable_h
#define DataTable_h 1

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>

class DataTableRowIterator;
class DataTableHeader;

class BadFormatException: public std::runtime_error
{
public:
	BadFormatException(const std::string& what_arg) :
		runtime_error(what_arg)
	{
	}
};

class WrongTypeException: public std::runtime_error
{
public:
	WrongTypeException(const std::string& what_arg) :
		runtime_error(what_arg)
	{
	}
};

/** @brief A data structure for holding named type columns
 *
 * DataTable is designed to hold tables of data where each column has a
 * name and can be of type integer, double or string.
 *
 * DataTable can be used to read and write TFS files, see DataTableReaderTFS
 * and DataTableWriterTFS.
 *
 * Note that accessing values requires typed versions of functions, as it
 * not possible to overload a function by return type. Use Get_d(), Get_i()
 * or Get_s(). For setting values the type can be inferred from the argument
 * type, Set(), or typed versions can be used Set_d(), Set_i() and Set_s().
 *
 * DataTable can also hold a set of values, also with types, in its header.
 */
class DataTable
{
protected:

	/// Used to index the location of columns in DataTable.
	struct location
	{
		char type;
		size_t pos;

	};

	//data storage
	std::vector<std::vector<double> > data_d;
	std::vector<std::vector<int> > data_i;
	std::vector<std::vector<std::string> > data_s;

	//indexing
	std::vector<std::string> col_names;
	std::unordered_map<std::string, location> lookup;
	size_t length;

	// Header values
	std::vector<double> hdata_d;
	std::vector<int> hdata_i;
	std::vector<std::string> hdata_s;
	std::vector<std::string> hcol_names;
	std::unordered_map<std::string, location> hlookup;

public:
	/// Construct an empty DataTable.
	DataTable() :
		length()
	{
	}

	/// Add a new column.
	/// \param type One of 'i' (integer), 'd' (double) or 's' (string)
	void AddColumn(std::string col_name, char type);
	/// New empty row.
	size_t AddRow();

	// get type in name
	///Get double value by name and row.
	double Get_d(const std::string col_name, size_t i) const;
	///Get integer value by name and row.
	int Get_i(const std::string col_name, size_t i) const;
	///Get string value by name and row.
	std::string Get_s(const std::string col_name, size_t i) const;
	///Get value by name and row converted to string.
	std::string GetAsStr(const std::string col_name, size_t i) const;

	// overloaded
	/// Set double value by name and row.
	void Set(const std::string col_name, size_t i, double x);
	/// Set integer value by name and row.
	void Set(const std::string col_name, size_t i, int x);
	/// Set string value by name and row.
	void Set(const std::string col_name, size_t i, std::string x);

	// set with type in name
	/// Set double value by name and row.
	void Set_d(const std::string col_name, size_t i, double x)
	{
		Set(col_name, i, x);
	}
	/// Set integer value by name and row.
	void Set_i(const std::string col_name, size_t i, int x)
	{
		Set(col_name, i, x);
	}
	/// Set string value by name and row.
	void Set_s(const std::string col_name, size_t i, std::string x)
	{
		Set(col_name, i, x);
	}
	/// Set value by name and row, converted from a string.
	void SetWithStr(const std::string col_name, size_t i, std::string x);

	/** @brief Variadic method for setting a whole row.
	 *
	 * Pass a set of values with the correct types for the columns, e.g.:
	 *
	 *     dt.AddRow('x', 1, 2.3, 4.2);
	 */
	template<typename ... Args>
	void AddRow(Args ... arg);

	DataTableRowIterator begin() const;
	DataTableRowIterator end() const;

	/// Access column names.
	const std::vector<std::string>& ColumnNames()
	{
		return col_names;
	}
	/// Access header names.
	const std::vector<std::string>& HeaderNames()
	{
		return hcol_names;
	}
	/// Access column types.
	char GetColumnType(std::string col_name)
	{
		return lookup.at(col_name).type;
	}
	/// Access header types.
	char GetHeaderType(std::string header_name)
	{
		return hlookup.at(header_name).type;
	}

	/// Add new element to the header .
	void HeaderAddColumn(std::string col_name, char type);

	// header access
	/// Get double value from header by name.
	double HeaderGet_d(const std::string col_name) const;
	/// Get integer value from header by name.
	int HeaderGet_i(const std::string col_name) const;
	/// Get string value from header by name.
	std::string HeaderGet_s(const std::string col_name) const;
	/// Get value from header by name as string.
	std::string HeaderGetAsStr(const std::string col_name) const;

	// overloaded
	/// Set double header value.
	void HeaderSet(const std::string col_name, double x);
	/// Set integer header value.
	void HeaderSet(const std::string col_name, int x);
	/// Set string header value.
	void HeaderSet(const std::string col_name, std::string x);

	// get type in name
	/// Set double header value.
	void HeaderSet_d(const std::string col_name, double x)
	{
		HeaderSet(col_name, x);
	}
	/// Set integer header value.
	void HeaderSet_i(const std::string col_name, int x)
	{
		HeaderSet(col_name, x);
	}
	/// Set string header value.
	void HeaderSet_s(const std::string col_name, std::string x)
	{
		HeaderSet(col_name, x);
	}
	/// Set value from header by name as string.
	void HeaderSetWithStr(const std::string col_name, std::string x);

	/// Return the length of table, i.e. number of rows.
	size_t Length() const
	{
		return length;
	}

	/// Test if the DataTable has a given column
	bool HasCol(std::string col_name) const;
	/// Test if the DataTable has a given column of given type
	bool HasCol(std::string col_name, char type) const;

	/// Test if the DataTable header has a given key
	bool HeaderHasKey(std::string col_name) const;
	/// Test if the DataTable header has a given key of given type
	bool HeaderHasKey(std::string col_name, char type) const;

private:
	/// See AddRow()
	template<typename T, typename ... Args>
	void AddRowN(size_t col_n, size_t row_n, T x, Args ... arg);
	/// See AddRow()
	void AddRowN(size_t, size_t)
	{
	}                               // Terminating case

public:
	//demo output
	/// Output data table to ascii stream.
	void OutputAscii(std::ostream &os) const;
};

template<typename ... Args>
void DataTable::AddRow(Args ... arg)
{
	size_t col_n = 0;
	size_t row_n = AddRow();
	AddRowN(col_n, row_n, arg ...);
}

template<typename T, typename ... Args>
void DataTable::AddRowN(size_t col_n, size_t row_n, T x, Args ... arg)
{
	std::string col_name = col_names.at(col_n);

	Set(col_name, row_n, x);
	AddRowN(col_n + 1, row_n, arg ...);
}

/** @brief Access individual rows of a DataTable.
 *
 * Type returned when using DataTableRowIterator. Does not actually hold
 * any data, just access values from the DataTable that is being iterated.
 *
 * Methods are wrappers around the accessors in DataTable, but without the
 * row number argument.
 */
class DataTableRow
{
public:
	DataTableRow(const DataTable *_dt, size_t _pos) :
		dt(_dt), pos(_pos)
	{
	}

	double Get_d(const std::string col_name) const
	{
		return dt->Get_d(col_name, pos);
	}
	int Get_i(const std::string col_name) const
	{
		return dt->Get_i(col_name, pos);
	}
	std::string Get_s(const std::string col_name) const
	{
		return dt->Get_s(col_name, pos);
	}
	std::string GetAsStr(const std::string col_name) const
	{
		return dt->GetAsStr(col_name, pos);
	}

private:
	const DataTable * dt;
	size_t pos;
};

/** @brief Emulate pointer access to individual rows of a DataTable
 *
 * See DataTableRow
 */
class DataTableRowPtr
{
public:
	DataTableRowPtr(const DataTable *_dt, size_t _pos) :
		dtr(_dt, _pos)
	{
	}
	DataTableRow &operator*()
	{
		return dtr;
	}
	DataTableRow *operator->()
	{
		return &dtr;
	}
private:
	DataTableRow dtr;
};

/// @brief Row iterator for DataTable.
class DataTableRowIterator: public std::iterator<std::bidirectional_iterator_tag,
	DataTableRow,
	std::ptrdiff_t,
	DataTableRow*,
	DataTableRow&>
{
public:
	DataTableRowIterator(const DataTable *_dt, size_t _pos) :
		dt(_dt), pos(_pos)
	{
	}
	DataTableRow operator *();
	DataTableRowPtr operator ->();
	DataTableRowIterator& operator++();
	DataTableRowIterator operator++(int);
	DataTableRowIterator& operator--();
	DataTableRowIterator operator--(int);
	bool operator==(const DataTableRowIterator &other) const;
	bool operator!=(const DataTableRowIterator &other) const;

private:
	const DataTable * dt;
	size_t pos;
};

#endif
