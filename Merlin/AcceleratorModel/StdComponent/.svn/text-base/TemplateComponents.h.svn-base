/*
 * Merlin C++ Class Library for Charged Particle Accelerator Simulations
 * 
 * Class library version 2.0 (2000)
 * 
 * file Merlin\AcceleratorModel\StdComponent\TemplateComponents.h
 * last modified 03/04/01 14:48:48
 */

/*
 * This file is derived from software bearing the following
 * restrictions:
 *
 * MERLIN C++ class library for 
 * Charge Particle Accelerator Simulations
 *
 * Copyright (c) 2000 by The Merlin Collaboration.  
 * ALL RIGHTS RESERVED. 
 *
 * Permission to use, copy, modify, distribute and sell this
 * software and its documentation for any purpose is hereby
 * granted without fee, provided that the above copyright notice
 * appear in all copies and that both that copyright notice and
 * this permission notice appear in supporting documentation.
 * No representations about the suitability of this software for
 * any purpose is made. It is provided "as is" without express
 * or implied warranty.
 */


#ifndef TemplateComponents_h
#define TemplateComponents_h 1

#include "merlin_config.h"


// AcceleratorComponent
#include "AcceleratorModel/AcceleratorComponent.h"



//	Template class for generating general accelerator
//	components with geometry type G and no field.

template <class G>
class TAccCompG : public AcceleratorComponent
{
public:

    typedef G geom_type;

public:
    TAccCompG (const string& id, G* geom, EMField* field = 0);

    TAccCompG (const TAccCompG<G>& rhs);


    TAccCompG<G>& operator = (const TAccCompG<G>& rhs);

    G& GetGeometry ();

    const G& GetGeometry () const;

protected:
private:
private:
};

//	Template class for generating general accelerator
//	components with geometry type G and field type F.

template <class G, class F>
class TAccCompGF : public TAccCompG<G>
{
public:

    typedef F field_type;

public:
    TAccCompGF (const string& id, G* geom, F* field);

    TAccCompGF (const TAccCompGF<G,F>& rhs);


    TAccCompGF<G,F>& operator = (const TAccCompGF<G,F>& rhs);

    F& GetField ();

    const F& GetField () const;

protected:
private:
private:
};

// This template class does provide a copy mechanism.
template <class G, class F>
class TAccCompGF_NC : public TAccCompG<G>
{
public:

    typedef F field_type;

public:
    TAccCompGF_NC (const string& id, G* geom, F* field);
	
	// Field accessors
    F& GetField ();
    const F& GetField () const;

protected:
private:
private:
};


// Parameterized Class TAccCompG

template <class G>
inline TAccCompG<G>::TAccCompG (const string& id, G* geom, EMField* field)
        : AcceleratorComponent(id,geom,field)
{
}



template <class G>
inline G& TAccCompG<G>::GetGeometry ()
{
    return static_cast<G&>(*itsGeometry);
}

template <class G>
inline const G& TAccCompG<G>::GetGeometry () const
{
    return static_cast<const G&>(*itsGeometry);
}

// Parameterized Class TAccCompGF

template <class G, class F>
inline TAccCompGF<G,F>::TAccCompGF (const string& id, G* geom, F* field)
        : TAccCompG<G>(id,geom,field)
{
}

template <class G, class F>
inline F& TAccCompGF<G,F>::GetField ()
{
    return static_cast<F&>(*this->itsField);
}

template <class G, class F>
inline const F& TAccCompGF<G,F>::GetField () const
{
    return static_cast<const F&>(*this->itsField);
}

// Parameterized Class TAccCompG

template <class G>
TAccCompG<G>::TAccCompG (const TAccCompG<G>& rhs)
        : AcceleratorComponent(rhs.GetName(),new G(rhs.GetGeometry()),0)
{
}



template <class G>
TAccCompG<G>& TAccCompG<G>::operator = (const TAccCompG<G>& rhs)
{
    if(this!=&rhs) {
        SetName(rhs.GetName());
        GetGeometry() = rhs.GetGeometry();
    }
    return *this;
}

// Parameterized Class TAccCompGF

template <class G, class F>
TAccCompGF<G,F>::TAccCompGF (const TAccCompGF<G,F>& rhs)
        : TAccCompG<G>(rhs.GetName(),new G(rhs.GetGeometry()), new F(rhs.GetField()))
{
}



template <class G, class F>
TAccCompGF<G,F>& TAccCompGF<G,F>::operator = (const TAccCompGF<G,F>& rhs)
{
    if(this!=&rhs) {
        TAccCompG<G>::operator=(rhs);
        GetField() = rhs.GetField();
    }
    return *this;
}

// Parameterized Class TAccCompGF_NC

template <class G, class F>
inline TAccCompGF_NC<G,F>::TAccCompGF_NC (const string& id, G* geom, F* field)
        : TAccCompG<G>(id,geom,field)
{}

template <class G, class F>
inline F& TAccCompGF_NC<G,F>::GetField ()
{
    return static_cast<F&>(*this->itsField);
}

template <class G, class F>
inline const F& TAccCompGF_NC<G,F>::GetField () const
{
    return static_cast<const F&>(*this->itsField);
}


#endif
