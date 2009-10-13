/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/03/20 13:42:54 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef ModelElement_h
#define ModelElement_h 1

#include "merlin_config.h"
#include <string>
#include <vector>

//	Root class for all elements/components which are used to
//	construct an AcceleratorModel. All ModelElement objects
//	are characterised by a type string and an identifier
//	(name). The type string identifies the type of element
//	(class), while the name is specific to an instance of
//	that element type.

class ModelElement
{
public:
    //	Constructor taking the name of the element.
	explicit ModelElement (const std::string& aName ="<UNNAMED>");

    virtual ~ModelElement ();

    //	Return the name of the element.
	virtual const std::string& GetName () const;

    //	Return the type string for the element.
	virtual const std::string& GetType () const = 0;

    //	Return the qualified name of the component. The
    //	qualified name has the form typestr.namestr.
	std::string GetQualifiedName () const;

    //	Virtual constructor.
    virtual ModelElement* Copy () const = 0;

    //	Set the name of the component.
	void SetName (const std::string& name);

	// Returns in ivec an ordered list of beamline indecies
	// associated with this ModelElement. Returns the length
	// of ivec.
	size_t GetBeamlineIndecies(std::vector<size_t>& ivec) const;
	virtual void AppendBeamlineIndecies(std::vector<size_t>& ivec) const =0;

protected:

    //	Initialise the ModelElement with the specified name.
	void Init (const std::string& aName);
	std::string id;
};

// Class ModelElement

inline ModelElement::ModelElement (const std::string& aName)
        : id(aName)
{
}

inline ModelElement::~ModelElement ()
{
    // nothing to do
}

inline const std::string& ModelElement::GetName () const
{
    return id;
}

inline std::string ModelElement::GetQualifiedName () const
{
    return GetType()+'.'+GetName();
}

inline void ModelElement::SetName (const std::string& name)
{
    id=name;
}

inline void ModelElement::Init (const std::string& aName)
{
    id=aName;
}

inline size_t ModelElement::GetBeamlineIndecies(std::vector<size_t>& ivec) const
{
	ivec.clear();
	AppendBeamlineIndecies(ivec);
	return ivec.size();
}

// utility macros:
// GetType() implementation
#define _TYPESTR(s) static const std::string typestr(#s); return typestr;

#endif
