/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "DataTable.h"
#include "../tests.h"
#include <iostream>
using namespace std;

/*
 * Test some of the basic features of DataTable
 *
 */

DataTable make_example_dt()
{
	DataTable dt1;
	dt1.AddColumn("a", 's');
	dt1.AddColumn("b", 'd');
	dt1.AddColumn("c", 'i');
	dt1.AddRow("alpha", 1.1, 1);
	dt1.AddRow("beta", 2.1, -2);
	dt1.AddRow("gamma", 3.1, 100);
	return dt1;
}

DataTable make_example_dt2()
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

void test1()
{
	cout << "test1()" << endl;
	DataTable dt1;
	dt1.AddColumn("name", 's');
	dt1.AddColumn("x", 'd');
	dt1.AddColumn("y", 'd');
	dt1.AddColumn("a", 'i');
	dt1.AddRow();

	assert(dt1.Get_d("x", 0) == 0);

	dt1.Set_d("x", 0, 5);
	assert(dt1.Get_d("x", 0) == 5);

	dt1.Set("name", 0, "Q1");
	assert(dt1.Get_s("name", 0) == "Q1");

	dt1.AddRow("D1", 1.1, 2.2, 5);
	assert(dt1.Get_s("name", 1) == "D1");
	assert(dt1.Get_d("x", 1) == 1.1);
	assert(dt1.Get_d("y", 1) == 2.2);
	assert(dt1.Get_i("a", 1) == 5);

	assert(dt1.GetAsStr("name", 1) == "D1");
	assert(dt1.GetAsStr("x", 1) == "1.100000");
	assert(dt1.GetAsStr("y", 1) == "2.200000");
	assert(dt1.GetAsStr("a", 1) == "5");

	assert(dt1.Get_d("x", 0) == 5);

	dt1.AddRow();
	dt1.SetWithStr("name", 2, "Q1");
	dt1.SetWithStr("x", 2, "7.7");
	dt1.SetWithStr("y", 2, "8.8");
	dt1.SetWithStr("a", 2, "9");
}

void test_get_wrong_col_type(DataTable &dt)
{
	cout << "test_get_wrong_col_type()" << endl;
	assert_throws(dt.Get_d("a", 0), WrongTypeException);
	assert_throws(dt.Get_d("a", 1), WrongTypeException);
	assert_throws(dt.Get_d("c", 0), WrongTypeException);

	assert_throws(dt.Get_i("a", 0), WrongTypeException);
	assert_throws(dt.Get_i("b", 0), WrongTypeException);

	assert_throws(dt.Get_s("b", 0), WrongTypeException);
	assert_throws(dt.Get_s("c", 0), WrongTypeException);

	assert_throws(dt.Set_d("a", 0, 0), WrongTypeException);
	assert_throws(dt.Set_d("a", 1, 0), WrongTypeException);
	assert_throws(dt.Set_d("c", 0, 0), WrongTypeException);
	assert_throws(dt.Set_i("b", 0, 0), WrongTypeException);

	assert_throws(dt.SetWithStr("b", 0, "a"), std::invalid_argument);

}

void test_range_errors(DataTable &dt)
{
	cout << "test_range_errors()" << endl;
	assert_throws(dt.Get_s("a", 3), std::out_of_range);
	assert_throws(dt.Get_d("b", 3), std::out_of_range);
	assert_throws(dt.Get_i("c", 3), std::out_of_range);
	assert_throws(dt.Get_i("d", 1), std::out_of_range);
}

void test_bad_add_col()
{
	cout << "test_bad_add_col()" << endl;
	DataTable dt1;
	dt1.AddColumn("a", 's');
	dt1.AddColumn("b", 'd');
	dt1.AddColumn("c", 'd');
	assert_throws(dt1.AddColumn("d", 'z'), WrongTypeException);
	assert_throws(dt1.AddColumn("c", 'd'), std::invalid_argument);
}

void test_add_col(DataTable dt)
{
	cout << "test_add_col()" << endl;

	dt.AddColumn("d", 'd');
	dt.AddColumn("e", 'i');
	dt.AddColumn("f", 's');
	assert(dt.Get_d("d", 0) == 0);
	assert(dt.Get_i("e", 0) == 0);
	assert(dt.Get_s("f", 0) == "");
}

void test_const(const DataTable dt)
{
	cout << "test_const()" << endl;
	assert(dt.Get_d("b", 0) == 1.1);
	dt.OutputAscii(std::cout);
	assert(dt.HeaderGet_s("z") == "test");
}

void test_has_key()
{
	cout << "test_has_key()" << endl;
	const auto dt2 = make_example_dt2();
	assert(dt2.HasCol("a") == true);
	assert(dt2.HasCol("d") == false);

	assert(dt2.HasCol("a", 's') == true);
	assert(dt2.HasCol("a", 'i') == false);
	assert(dt2.HasCol("d", 'i') == false);

	assert(dt2.HeaderHasKey("x") == true);
	assert(dt2.HeaderHasKey("w") == false);

	assert(dt2.HeaderHasKey("x", 'd') == true);
	assert(dt2.HeaderHasKey("x", 'i') == false);
	assert(dt2.HeaderHasKey("w", 'i') == false);
}

int main()
{
	test1();

	auto dt1 = make_example_dt();
	test_get_wrong_col_type(dt1);
	test_range_errors(dt1);
	test_bad_add_col();
	test_add_col(dt1);

	const auto dt2 = make_example_dt2();
	test_const(dt2);

	test_has_key();

	return 0;
}
