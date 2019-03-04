/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef AMBufferManager_h
#define AMBufferManager_h 1

#include "merlin_config.h"
#include <set>
#include <algorithm>
#include <cassert>

/*
 *	Template class for managing active diagnostic buffers.
 *	The template M should be a diagnostic type defining a
 *	type M::Data and a type M::Buffer. M::Buffer should
 *	supply the following method:
 */

//	void M::Buffer::Record(const M&, const M::Data&)

template<class M, class B, class D>
class AMBufferManager
{
public:

	/*
	 *	Add the specified buffer.
	 */
	void AddBuffer(B* buf);

	/*
	 *	Remove buf from the buffer list, if it exists. Returns
	 *	true if successful.
	 */
	bool RemoveBuffer(B* buf);

	/*
	 *	Remove all buffers (not including the default buffer).
	 */
	void ClearAllBuffers();

	/*
	 *	Sets the default buffer for all diagnostics of type M.
	 */
	static void SetDefaultBuffer(B* buf);

	/*
	 *	Sends the data to all the buffers.
	 */
	void SendToBuffers(const M& monitor, const D& data);

	/*
	 *	Returns true if there are no buffers.
	 *	@return True if there are no buffers
	 */
	bool empty() const;

private:

	std::set<B*> buffers;
	static B* defBuffer;
};

template<class M, class B, class D>
inline bool AMBufferManager<M, B, D>::empty() const
{
	return defBuffer == nullptr && buffers.empty();
}

template<class M, class B, class D>
B * AMBufferManager<M, B, D>::defBuffer = nullptr;

template<class M, class B, class D>
void AMBufferManager<M, B, D>::AddBuffer(B* buf)
{
	buffers.insert(buf);
}

template<class M, class B, class D>
bool AMBufferManager<M, B, D>::RemoveBuffer(B* buf)
{
#ifndef NDEBUG
	int n = buffers.erase(buf);
	assert(n == 0 || n == 1);
	return n == 1;
#else
	return buffers.erase(buf) != 0;
#endif
}

template<class M, class B, class D>
void AMBufferManager<M, B, D>::ClearAllBuffers()
{
	buffers.clear();
}

template<class M, class B, class D>
void AMBufferManager<M, B, D>::SetDefaultBuffer(B* buf)
{
	defBuffer = buf;
}

template<class M, class B, class D>
void AMBufferManager<M, B, D>::SendToBuffers(const M& monitor, const D& data)
{
	if(defBuffer != nullptr)
	{
		defBuffer->Record(monitor, data);
	}
	for(typename std::set<B*>::iterator b = buffers.begin(); b != buffers.end(); b++)
	{
		(*b)->Record(monitor, data);
	}
}

#endif
