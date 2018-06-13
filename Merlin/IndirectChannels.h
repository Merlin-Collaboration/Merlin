/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef IndirectChannels_h
#define IndirectChannels_h 1

#include "merlin_config.h"
#include "ChannelServer.h"
#include "Channels.h"

template<class E> class TIRWChannel;

/**
 *	Template class for constructing TIRWChannel objects for
 *	specific ModelElement types.
 */

template<class E>
class TIC_ctor: public ChannelServer::ChannelCtor
{
public:

	typedef E element_type;
	typedef double (E::* read_f)() const;
	typedef void (E::* write_f)(double);

	TIC_ctor(const string& aType, const string& aKey, read_f rmf, write_f wmf = 0);

	/**
	 *	Constructs a channel for the specified ModelElement.
	 */
	virtual ROChannel* ConstructRO(ModelElement* anElement);

	/**
	 *	Constructs a channel for the specified ModelElement.
	 */
	virtual RWChannel* ConstructRW(ModelElement* anElement);

	double ReadFrom(const E* elmnt);
	void WriteTo(E* elmnt, double value);

private:

	read_f r_f;
	write_f w_f;
};

/**
 *	Template class for constructing a concrete Indirect
 *	RWChannel object. E can be any concrete ModelElement
 *	type, which has set- and get-like member functions to
 *	the associated parameter. The functions must have the
 *	following signature:
 *
 *	get: double (E::*)() const
 *	set: void (E::*)(double)
 */

template<class E>
class TIRWChannel: public RWChannel
{
public:

	typedef TIC_ctor<E> ch_ctor;

	TIRWChannel(E* element, TIC_ctor<E>* proto_f);

	/**
	 *	Returns the ID of the channel (parameter).
	 *	@return Channel ID
	 */
	virtual string GetID() const;

	/**
	 *	Returns the current value of the parameter/attribute
	 *	associated with the channel.
	 *
	 *	@return Value of the channel's associated parameter/attribute
	 */
	virtual double Read() const;

	/**
	 *	Sets the current value of the parameter/attribute
	 *	associated with the channel.
	 */
	virtual void Write(double value);

	/**
	 *	Increments the current value of the parameter/attribute
	 *	associated with the channel. Returns the final
	 *	(incremented) value.
	 *
	 *	@return Final (incremented) value
	 */
	virtual double Increment(double delta);

private:

	E* theElement;
	ch_ctor* fp;

	TIRWChannel(const TIRWChannel& rhs);
	TIRWChannel& operator=(const TIRWChannel& rhs);

};

template<class E>
TIC_ctor<E>::TIC_ctor(const string& aType, const string& aKey, read_f rmf, write_f wmf) :
	ChannelServer::ChannelCtor(aType, aKey), r_f(rmf), w_f(wmf)
{
}

template<class E>
ROChannel* TIC_ctor<E>::ConstructRO(ModelElement* anElement)
{
	return ConstructRW(anElement);
}

template<class E>
RWChannel* TIC_ctor<E>::ConstructRW(ModelElement* anElement)
{
#ifndef NDEBUG
	E* eptr = dynamic_cast<E*>(anElement);
	assert(eptr); // bad cast!
	return new TIRWChannel<E>(static_cast<E*>(eptr), this);
#else
	return new TIRWChannel<E>(static_cast<E*>(anElement), this);
#endif
}

template<class E>
double TIC_ctor<E>::ReadFrom(const E* elmnt)
{
	return (elmnt->*r_f)();
}

template<class E>
void TIC_ctor<E>::WriteTo(E* elmnt, double value)
{
	(elmnt->*w_f)(value);
}

template<class E>
TIRWChannel<E>::TIRWChannel(E* element, TIC_ctor<E>* proto_f) :
	theElement(element), fp(proto_f)
{
}

template<class E>
string TIRWChannel<E>::GetID() const
{
	return theElement->GetQualifiedName() + '.' + fp->GetKey();
}

template<class E>
double TIRWChannel<E>::Read() const
{
	return fp->ReadFrom(theElement);
}

template<class E>
void TIRWChannel<E>::Write(double value)
{
	fp->WriteTo(theElement, value);
}

template<class E>
double TIRWChannel<E>::Increment(double delta)
{
	double v = fp->ReadFrom(theElement) + delta;
	fp->WriteTo(theElement, v);
	return v;
}

#endif
