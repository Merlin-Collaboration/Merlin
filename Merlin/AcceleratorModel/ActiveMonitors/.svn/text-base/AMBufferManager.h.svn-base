/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/07 09:14:12 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef AMBufferManager_h
#define AMBufferManager_h 1

#include "merlin_config.h"
#include <set>
#include <algorithm>
#include <cassert>

//	Template class for managing active diagnotic buffers.
//	The template M should be a diagnostic type defining a
//	type M::Data and a type M::Buffer. M::Buffer should
//	supply the following method:
//
//	void M::Buffer::Record(const M&, const M::Data&)

template <class M, class B, class D>
class AMBufferManager
{
public:

    //	Add the specified buffer.
    void AddBuffer (B* buf);

    //	Remove buf from the buffer list, if it exisits. Returns
    //	true if successful.
    bool RemoveBuffer (B* buf);

    //	Remove all buffers (not including the default buffer).
    void ClearAllBuffers ();

    //	Sets the default buffer for all diagnostics of type M.
    static void SetDefaultBuffer (B* buf);

    //	Sends the data to all the buffers.
    void SendToBuffers (const M& monitor, const D& data);

    //	Returns true if there are no buffers.
    bool empty () const;

private:

    std::set<B*> buffers;
    static B* defBuffer;
};

template <class M, class B, class D>
inline bool AMBufferManager<M,B,D>::empty () const
{
    return defBuffer==0 && buffers.empty();
}

template <class M, class B, class D>
B* AMBufferManager<M,B,D>::defBuffer = 0;


template <class M, class B, class D>
void AMBufferManager<M,B,D>::AddBuffer (B* buf)
{
    buffers.insert(buf);
}

template <class M, class B, class D>
bool AMBufferManager<M,B,D>::RemoveBuffer (B* buf)
{
#ifndef NDEBUG
    int n = buffers.erase(buf);
    assert(n==0||n==1);
    return n==1;
#else
    return buffers.erase(buf)!=0;
#endif
}

template <class M, class B, class D>
void AMBufferManager<M,B,D>::ClearAllBuffers ()
{
    buffers.clear();
}

template <class M, class B, class D>
void AMBufferManager<M,B,D>::SetDefaultBuffer (B* buf)
{
    defBuffer = buf;
}

template <class M, class B, class D>
void AMBufferManager<M,B,D>::SendToBuffers (const M& monitor, const D& data)
{
    if(defBuffer!=0)
        defBuffer->Record(monitor,data);
    for(typename std::set<B*>::iterator b=buffers.begin(); b!=buffers.end(); b++)
            (*b)->Record(monitor,data);
}

#endif
