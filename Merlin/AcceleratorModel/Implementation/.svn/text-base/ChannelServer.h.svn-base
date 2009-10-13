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
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef ChannelServer_h
#define ChannelServer_h 1

#include "merlin_config.h"
#include <string>
#include <map>
#include "AcceleratorModel/AcceleratorModel.h"
#include "AcceleratorModel/Implementation/ElementRepository.h"

class ModelElement;
class StringPattern;
class RWChannel;
class ROChannel;

//	Responsible for providing an interface to ROChannel and
//	RWChannel for the ModelElements in the AcceleratorModel.

class ChannelServer
{
public:
    //	Abstract factory class for RWChannel objects.
    class ChannelCtor
    {
    public:
        // destructor
        virtual ~ChannelCtor() {}

        //	Constructs a channel for the specified ModelElement.
        virtual ROChannel* ConstructRO (ModelElement* anElement) = 0;

        //	Constructs a channel for the specified ModelElement.
        virtual RWChannel* ConstructRW (ModelElement* anElement) = 0;

        //	Returns the ID of the channel (i.e. type.key).
        std::string GetID ();

        //	Returns the ModelElement type.
        const string& GetType () const;

        //	Returns the channel key.
        const string& GetKey () const;

    protected:

        ChannelCtor (const string& aType, const string& aKey);

        std::string type;
        std::string key;
    };

    typedef std::map<std::string,ChannelCtor*> CtorTable;

    //	Returns in channels all ROChannels matching chID.
    //	Returns the number of channels found.
    size_t GetROChannels (const string& chID, std::vector<ROChannel*>& channels);

    //	Returns in channels all RWChannels matching chID.
    //	Returns the number of channels found.
    size_t GetRWChannels (const string& chID, std::vector<RWChannel*>& channels);

    //	Returns read-only channels matching chid for all
    //	matching components in aBeamline. Note that only
    //	channels associated with AcceleratorComponents can be
    //	extracted using this method.
    size_t GetROChannels (AcceleratorModel::Beamline& aBeamline, const std::string& chid, std::vector<ROChannel*>& channels);

    //	Returns read-write channels matching chid for all
    //	matching components in aBeamline. Note that only
    //	channels associated with AcceleratorComponents can be
    //	extracted using this method.
    size_t GetRWChannels (AcceleratorModel::Beamline& aBeamline, const std::string& chid, std::vector<RWChannel*>& channels);

    //	Adds a ChannelCtor object to the server.
    void RegisterCtor (ChannelCtor* chctor);

    void SetRepository (ElementRepository* me_repo);

private:

    ElementRepository* theElements;
    CtorTable chCtors;

    void FindCtors (const string& type, const string& keypat, std::set<ChannelCtor*>& ctors);

    //	Returns in elements the ModelElements that match
    //	pattern. elements is sorted by TYPE.
    void FindElements (const std::string& id_pat, std::vector<ModelElement*>& elements);
};

inline ChannelServer::ChannelCtor::ChannelCtor (const string& aType, const string& aKey)
        : type(aType),key(aKey)
{}

inline const string& ChannelServer::ChannelCtor::GetType () const
{
    return type;
}

inline const string& ChannelServer::ChannelCtor::GetKey () const
{
    return key;
}

#endif
