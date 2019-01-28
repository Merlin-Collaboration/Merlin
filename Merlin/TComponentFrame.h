/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef TComponentFrame_h
#define TComponentFrame_h 1

#include "merlin_config.h"
#include <iostream>
#include "ComponentFrame.h"

/**
 *	A Typed ComponentFrame. Used to construct Accelerator
 *	Component type specific ComponentFrame objects.
 */
template<class T>
class TComponentFrame: public ComponentFrame
{
public:

	typedef T ComponentType;

	/**
	 *	Constructor.
	 */
	explicit TComponentFrame(T& component);

	T& GetComponent();
	const T& GetComponent() const;
	T* operator ->();
	const T* operator ->() const;

	/**
	 *	Returns a copy of this ComponentFrame. Note that only
	 *	the reference to the AcceleratorComponent is copied, not
	 *	the AcceleratorComponent itself.
	 *
	 *	@return Copy of the ComponentFrame
	 */
	virtual ModelElement* Copy() const;
};

template<class T>
inline TComponentFrame<T>::TComponentFrame(T& component) :
	ComponentFrame(component)
{
}

template<class T>
inline T& TComponentFrame<T>::GetComponent()
{
	return *static_cast<T*>(theComponent);
}

template<class T>
inline const T& TComponentFrame<T>::GetComponent() const
{
	return *static_cast<const T*>(theComponent);
}

template<class T>
inline T* TComponentFrame<T>::operator ->()
{
	return static_cast<T*>(theComponent);
}

template<class T>
inline const T* TComponentFrame<T>::operator ->() const
{
	return static_cast<const T*>(theComponent);
}

template<class T>
inline ModelElement* TComponentFrame<T>::Copy() const
{
	return new TComponentFrame<T>(*this);
}

#endif
