/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_ILCDFS_IO
#define _h_ILCDFS_IO 1

#include <iostream>
#include <string>

class dfs_trace
{

	static std::ostream* os;
	const bool do_output;

public:

	enum trace_level
	{
		error=0,
		warning=1,
		level_1=2,
		level_2=3,
		level_3=4

	};
	static trace_level verbosity;

	dfs_trace(trace_level l = level_3) :
		do_output((l <= verbosity) && (os != nullptr))
	{
	}

	static void set_trace_stream(std::ostream& anOs)
	{
		os = &anOs;
	}

	template<class T>
	friend dfs_trace operator<<(dfs_trace dfst, const T& t)
	{
		if(dfst.do_output)
		{
			(*(dfst.os)) << t;
		}
		return dfst;
	}
	// support for manipulators (e.g. endl)
	friend dfs_trace operator<<(dfs_trace dfst, std::ostream& (*mp)(std::ostream&))
	{
		if(dfst.do_output)
		{
			mp(*(dfst.os));
		}
		return dfst;
	}

};

#endif
