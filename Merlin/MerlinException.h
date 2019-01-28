/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef MerlinException_h
#define MerlinException_h 1

#include <string>

/**
 * Root class for all Merlin exceptions.
 */
class MerlinException
{
public:

	/**
	 * Constructor: Builds the MerlinException and sets the exception message.
	 * @param[in] s The exception message.
	 */
	explicit MerlinException(const std::string& s);

	/**
	 * Constructor: Builds the MerlinException.
	 */
	MerlinException();

	/**
	 * Virtual destructor.
	 */
	virtual ~MerlinException()
	{
	}

	/**
	 * Gets the exception message.
	 * @return The current exception message.
	 */
	const char* Msg() const;

protected:

	/**
	 * Sets the exception message.
	 * @see AppendMsg
	 * @param[in] s The exception message.
	 */
	void SetMsg(const std::string& s);

	/**
	 * Appends a string to the exception message.
	 * @see SetMsg
	 * @param[in] s The string to append to the exception message.
	 */
	void AppendMsg(const std::string& s);

private:

	/**
	 * String storage containing the exception message.
	 */
	std::string msg;
};

inline MerlinException::MerlinException(const std::string& s) :
	msg(s)
{
}

inline MerlinException::MerlinException() :
	msg()
{
}

inline const char* MerlinException::Msg() const
{
	return msg.c_str();
}

inline void MerlinException::SetMsg(const std::string& s)
{
	msg = s;
}

inline void MerlinException::AppendMsg(const std::string& s)
{
	msg += s;
}

#endif
