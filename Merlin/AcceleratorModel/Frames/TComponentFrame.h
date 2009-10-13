/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef TComponentFrame_h
#define TComponentFrame_h 1

#include "merlin_config.h"
#include <iostream>
// ComponentFrame
#include "AcceleratorModel/Frames/ComponentFrame.h"

//	A Typed ComponentFrame. Used to construct Accelerator
//	Component type specific ComponentFrame objects.

template <class T>
class TComponentFrame : public ComponentFrame
{
public:

    typedef T ComponentType;

    //	Constructor.
    explicit TComponentFrame (T& component);

    T& GetComponent ();
    const T& GetComponent () const;
    T* operator -> ();
    const T* operator -> () const;

    //	Returns a copy of this ComponentFrame. Note that only
    //	the reference to the AcceleratorComponent is copied, not
    //	the AcceleratorComponent itself.
    virtual ModelElement* Copy () const;
};

template <class T>
inline TComponentFrame<T>::TComponentFrame (T& component)
        : ComponentFrame(component)
{}

template <class T>
inline T& TComponentFrame<T>::GetComponent ()
{
    return *static_cast<T*>(theComponent);
}

template <class T>
inline const T& TComponentFrame<T>::GetComponent () const
{
    return *static_cast<const T*>(theComponent);
}

template <class T>
inline T* TComponentFrame<T>::operator -> ()
{
    return static_cast<T*>(theComponent);
}

template <class T>
inline const T* TComponentFrame<T>::operator -> () const
{
    return static_cast<const T*>(theComponent);
}

template <class T>
inline ModelElement* TComponentFrame<T>::Copy () const
{
    return new TComponentFrame<T>(*this);
}

#endif
