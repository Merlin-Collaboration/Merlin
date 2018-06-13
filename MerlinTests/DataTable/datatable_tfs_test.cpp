/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "DataTable.h"
#include "DataTableTFS.h"
#include "../tests.h"
#include <iostream>
#include <sstream>
using namespace std;

/*
 * Test reading and writing of DataTable TFS files
 *
 */

DataTable make_example_dt()
{
	DataTable dt;
	dt.AddColumn("a", 's');
	dt.AddColumn("b", 'd');
	dt.AddColumn("c", 'i');
	dt.AddRow("alpha", 1.1, 1);
	dt.AddRow("beta", 2.1, -2);
	dt.AddRow("gamma", 3.1, 100);
	dt.HeaderAddColumn("x", 'd');
	dt.HeaderAddColumn("y", 'i');
	dt.HeaderAddColumn("z", 's');
	dt.HeaderSet("x", 9.9);
	dt.HeaderSet("y", 9);
	dt.HeaderSet("z", "test");
	return dt;
}

void write_read(DataTable dt)
{
	cout << "write_read()" << endl;
	stringstream ss;

	DataTableWriterTFS(&ss).Write(dt);
	unique_ptr<DataTable> dt2(DataTableReaderTFS(&ss).Read());

	assert(dt.HeaderGet_d("x") == dt2->HeaderGet_d("x"));
	assert(dt.HeaderGet_i("y") == dt2->HeaderGet_i("y"));
	assert(dt.HeaderGet_s("z") == dt2->HeaderGet_s("z"));

	assert(dt.Length() == dt2->Length());
	for(size_t i = 0; i < dt.Length(); i++)
	{
		assert(dt.Get_s("a", i) == dt2->Get_s("a", i));
		assert(dt.Get_d("b", i) == dt2->Get_d("b", i));
		assert(dt.Get_i("c", i) == dt2->Get_i("c", i));
	}
}

void read_big()
{
	cout << "read_big()" << endl;
	string paths[] = {"../data/twiss.7.0tev.b1_new.tfs", "data/twiss.7.0tev.b1_new.tfs", "MerlinTests/data/twiss.7.0tev.b1_new.tfs"};

	string lattice_path;
	for(size_t i = 0; i < 3; i++)
	{
		ifstream test_file;
		test_file.open(paths[i].c_str());
		if(test_file)
		{
			lattice_path = paths[i];
			break;
		}
	}
	cout << "Reading: " << lattice_path << endl;

	unique_ptr<DataTable> dt(DataTableReaderTFS(lattice_path).Read());

	assert(dt->Length() == 13211);

	assert(dt->HeaderGet_s("NAME") == "TWISS");
	assert_close(dt->HeaderGet_d("MASS"), 9.3827201299999996E-01, 1e-8);
	assert(dt->Get_s("NAME", 0) == "LHCB1$START");
	assert(dt->Get_s("NAME", 2) == "MBAS2.1R1");
	assert_close(dt->Get_d("S", 4), 2.0915000000000003E+01, 1e-8);
}

int main()
{
	auto dt1 = make_example_dt();

	write_read(dt1);
	read_big();

	return 0;
}
