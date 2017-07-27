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
	BadFormatException(const std::string& what_arg): runtime_error(what_arg) {};
};

class WrongTypeException: public std::runtime_error
{
public:
	WrongTypeException(const std::string& what_arg): runtime_error(what_arg) {};
};

class DataTable
{
protected:
	struct location
	{
		char type;
		size_t pos;
	};

	//data storage
	std::vector<std::vector<double>> data_d;
	std::vector<std::vector<int>> data_i;
	std::vector<std::vector<std::string>> data_s;

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
	DataTable():length() {};

	void AddColumn(std::string col_name, char type);
	size_t AddRow();

	// get type in name
	double Get_d(const std::string col_name, size_t i) const;
	int Get_i(const std::string col_name, size_t i) const;
	std::string Get_s(const std::string col_name, size_t i) const;
	std::string GetAsStr(const std::string col_name, size_t i) const;

	// overloaded
	void Set(const std::string col_name, size_t i, double x);
	void Set(const std::string col_name, size_t i, int x);
	void Set(const std::string col_name, size_t i, std::string x);

	// get type in name
	void Set_d(const std::string col_name, size_t i, double x)
	{
		Set(col_name, i, x);
	}
	void Set_i(const std::string col_name, size_t i, int x)
	{
		Set(col_name, i, x);
	}
	void Set_s(const std::string col_name, size_t i, std::string x)
	{
		Set(col_name, i, x);
	}
	void SetWithStr(const std::string col_name, size_t i, std::string x);

	// add row variable length argument list
	template <typename... Args>
	void AddRow(Args... arg);

	DataTableRowIterator begin() const;
	DataTableRowIterator end() const;

	const std::vector<std::string>& ColumnNames()
	{
		return col_names;
	}
	const std::vector<std::string>& HeaderNames()
	{
		return hcol_names;
	}
	char GetColumnType(std::string col_name)
	{
		return lookup.at(col_name).type;
	}
	char GetHeaderType(std::string header_name)
	{
		return hlookup.at(header_name).type;
	}

	void HeaderAddColumn(std::string col_name, char type);

	// header access
	double HeaderGet_d(const std::string col_name) const;
	int HeaderGet_i(const std::string col_name) const;
	std::string HeaderGet_s(const std::string col_name) const;
	std::string HeaderGetAsStr(const std::string col_name) const;

	// overloaded
	void HeaderSet(const std::string col_name, double x);
	void HeaderSet(const std::string col_name, int x);
	void HeaderSet(const std::string col_name, std::string x);

	// get type in name
	void HeaderSet_d(const std::string col_name, double x)
	{
		HeaderSet(col_name, x);
	}
	void HeaderSet_i(const std::string col_name, int x)
	{
		HeaderSet(col_name, x);
	}
	void HeaderSet_s(const std::string col_name, std::string x)
	{
		HeaderSet(col_name, x);
	}
	void HeaderSetWithStr(const std::string col_name, std::string x);

	size_t Length() const
	{
		return length;
	};
private:
	template <typename T, typename... Args>
	void AddRowN(size_t col_n, size_t row_n, T x, Args... arg);
	void AddRowN(size_t, size_t) {}; // Terminating case

public:
	//demo output
	void OutputAscii(std::ostream &os) const;
};


template <typename... Args>
void DataTable::AddRow(Args... arg)
{
	size_t col_n = 0;
	size_t row_n = AddRow();
	AddRowN(col_n, row_n, arg...);
}

template <typename T, typename... Args>
void DataTable::AddRowN(size_t col_n, size_t row_n, T x, Args... arg)
{
	std::string col_name = col_names.at(col_n);

	Set(col_name, row_n, x);
	AddRowN(col_n+1, row_n, arg...);
}

class DataTableRow
{
public:
	DataTableRow(const DataTable *_dt, size_t _pos): dt(_dt), pos(_pos) {}

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

class DataTableRowPtr
{
public:
	DataTableRowPtr(const DataTable *_dt, size_t _pos): dtr(_dt, _pos) {}
	DataTableRow &operator* ()
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


class DataTableRowIterator :public std::iterator<std::bidirectional_iterator_tag,
	DataTableRow,
	std::ptrdiff_t,
	DataTableRow*,
	DataTableRow&>
{
public:
	DataTableRowIterator(const DataTable *_dt, size_t _pos): dt(_dt), pos(_pos) {}
	DataTableRow operator * ();
	DataTableRowPtr operator -> ();
	DataTableRowIterator& operator++ ();
	DataTableRowIterator operator++ (int);
	DataTableRowIterator& operator-- ();
	DataTableRowIterator operator-- (int);
	bool operator== (const DataTableRowIterator &other) const;
	bool operator!= (const DataTableRowIterator &other) const;

private:
	const DataTable * dt;
	size_t pos;
};

#endif
