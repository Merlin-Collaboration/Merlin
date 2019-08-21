/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef MerlinException_h
#define MerlinException_h 1

#include <string>
#include <stdexcept>

/**
 * Root class for all Merlin exceptions.
 */
class MerlinException: public std::runtime_error
{
public:
	/**
	 * Constructor: Builds the MerlinException and sets the exception message.
	 * @param[in] what_arg The exception message.
	 */
	MerlinException(const std::string& what_arg) :
		std::runtime_error(what_arg)
	{
	}

	/// Deprecated use MerlinException::what()
	[[deprecated("Use what()")]]
	const char* Msg() const noexcept
	{
		return what();
	}

	/**
	 * \fn virtual const char* MerlinException::what() const noexcept;
	 * \memberof MerlinException
	 *
	 * Returns the explanatory string. Inheritied from std::runtime_error
	 */
};

#endif
