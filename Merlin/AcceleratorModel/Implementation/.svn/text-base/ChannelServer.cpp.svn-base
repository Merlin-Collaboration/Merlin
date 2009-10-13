/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/04/26 20:02:47 $
// $Revision: 1.5 $
// 
/////////////////////////////////////////////////////////////////////////

#include <cassert>
#include "IO/MerlinIO.h"
// Channels
#include "Channels/Channels.h"
// ModelElement
#include "AcceleratorModel/ModelElement.h"
// StringPattern
#include "utility/StringPattern.h"
// ChannelServer
#include "AcceleratorModel/Implementation/ChannelServer.h"

namespace {
    
    struct CHPath {
        
        string type;
        string id;
        string key;
        
        CHPath(const string& id);
        
        string typekey() const {
            return type+'.'+key;
        }
        string eid() const {
            return type+'.'+id;
        }
    };
    
    CHPath::CHPath(const string& cid)
    {
        int n1=cid.find('.');
        type=cid.substr(0,n1);
        
        int n2=cid.find('.',n1+1);
        id=cid.substr(n1+1,n2-n1-1);
        
        key=cid.substr(n2+1);
    }
    
    
    struct TYPE_CMP {
        bool operator()(const ModelElement* e1, const ModelElement* e2) {
            return (*e1).GetType() < (*e2).GetType();
        }
    };
    
    struct FindType {
        string type;
        bool neg;
        
        FindType(const string& aType, bool negate=false)
            : type(aType),neg(negate) {}
        
        bool operator()(const ChannelServer::CtorTable::value_type& ctor) {
            bool ok = (ctor.second)->GetType() == type;
            return neg ? !ok : ok;
        }
    };
    
}; // end anonymous namespace

std::string ChannelServer::ChannelCtor::GetID ()
{
    return type+'.'+key;
}

size_t ChannelServer::GetROChannels (const string& chID, std::vector<ROChannel*>& channels)
{
    CHPath chpath(chID);
    vector<ModelElement*> elements;
    elements.reserve(8);
    
    FindElements(chpath.eid(),elements);
    if(elements.empty()) {
        MerlinIO::warning()<<"ChannelServer: No elements matching "<<chpath.eid()<<" found"<<endl;
        return 0;
    }
    
    set<ChannelCtor*> ctors;
    string etype="";
    int count=0;
    for(vector<ModelElement*>::iterator ei = elements.begin(); ei!=elements.end(); ei++) {
        
        // If we have a new type, get the corresponding constructors
        if(etype!=(*ei)->GetType()) {
            etype=(*ei)->GetType();
            FindCtors(etype,chpath.key,ctors);
            if(ctors.empty()) {
                MerlinIO::warning()<<"ChannelServer: No channel constructors for "<<chpath.type<<" found"<<endl;
                break;
            }
        }
        
        // Now construct the channels
        for(set<ChannelCtor*>::iterator ci = ctors.begin(); ci!=ctors.end(); ci++) {
            channels.push_back((*ci)->ConstructRO(*ei));
            count++;
        }
    }
    
    return count;
}

size_t ChannelServer::GetRWChannels (const string& chID, std::vector<RWChannel*>& channels)
{
    CHPath chpath(chID);
    vector<ModelElement*> elements;
    elements.reserve(8);
    
    FindElements(chpath.eid(),elements);
    if(elements.empty()) {
        MerlinIO::warning()<<"ChannelServer: No elements matching "<<chpath.eid()<<" found"<<endl;
        return 0;
    }
    
    set<ChannelCtor*> ctors;
    string etype="";
    int count=0;
    for(vector<ModelElement*>::iterator ei = elements.begin(); ei!=elements.end(); ei++) {
        
        // If we have a new type, get the corresponding constructors
        if(etype!=(*ei)->GetType()) {
            etype=(*ei)->GetType();
            FindCtors(etype,chpath.key,ctors);
            if(ctors.empty()) {
                MerlinIO::warning()<<"ChannelServer: No channel constructors for "<<chpath.type<<" found"<<endl;
                break;
            }
        }
        
        // Now construct the channels
        for(set<ChannelCtor*>::iterator ci = ctors.begin(); ci!=ctors.end(); ci++) {
            RWChannel* ch = (*ci)->ConstructRW(*ei);
            assert(ch!=0); // if(ch==0) throw ChannelIsRO(ch->GetID());
            channels.push_back(ch);
            count++;
        }
    }
    
    return count;
}

size_t ChannelServer::GetROChannels (AcceleratorModel::Beamline& aBeamline, const std::string& chid, std::vector<ROChannel*>& channels)
{
    CHPath chpath(chid);
    StringPattern idpat(chpath.eid());
    size_t count = 0;
    
    for(AcceleratorModel::BeamlineIterator cmpi = aBeamline.begin(); cmpi!=aBeamline.end(); cmpi++) {
        
        if((*cmpi)->IsComponent()) {
            AcceleratorComponent& component = (*cmpi)->GetComponent();
            
            if(idpat(component.GetQualifiedName())) {
                
                set<ChannelCtor*> ctors;
                FindCtors(component.GetType(),chpath.key,ctors);
                
                // Now construct the channels for the current component
                for(set<ChannelCtor*>::iterator ci = ctors.begin(); ci!=ctors.end(); ci++) {
                    ROChannel* ch = (*ci)->ConstructRO(&component);
                    assert(ch!=0);
                    channels.push_back(ch);
                    count++;
                }
            }
        }
    }
    
    return count;
}

size_t ChannelServer::GetRWChannels (AcceleratorModel::Beamline& aBeamline, const std::string& chid, std::vector<RWChannel*>& channels)
{
    CHPath chpath(chid);
    StringPattern idpat(chpath.eid());
    size_t count = 0;
    
    for(AcceleratorModel::BeamlineIterator cmpi = aBeamline.begin(); cmpi!=aBeamline.end(); cmpi++) {
        if((*cmpi)->IsComponent()) {
            
            AcceleratorComponent& component = (*cmpi)->GetComponent();
            
            if(idpat(component.GetQualifiedName())) {
                
                set<ChannelCtor*> ctors;
                FindCtors(component.GetType(),chpath.key,ctors);
                
                // Now construct the channels for the current component
                for(set<ChannelCtor*>::iterator ci = ctors.begin(); ci!=ctors.end(); ci++) {
                    RWChannel* ch = (*ci)->ConstructRW(&component);
                    assert(ch!=0);
                    channels.push_back(ch);
                    count++;
                }
            }
        }
    }
    
    return count;
}

void ChannelServer::RegisterCtor (ChannelCtor* chctor)
{
    pair<CtorTable::iterator,bool> rv = chCtors.insert(
        CtorTable::value_type(chctor->GetID(),chctor));
    assert(rv.second);
}

void ChannelServer::SetRepository (ElementRepository* me_repo)
{
    theElements = me_repo;
}

void ChannelServer::FindCtors (const string& type, const string& keypat, std::set<ChannelCtor*>& ctors)
{
    ctors.clear();
    
    CtorTable::iterator c1 = find_if(chCtors.begin(),chCtors.end(),FindType(type));
    
    if(c1==chCtors.end())
        return;
    
    CtorTable::iterator c2 = find_if(c1,chCtors.end(),FindType(type,true));
    
    // Now find those ctors that match the key pattern
    StringPattern keyp(keypat);
    for(CtorTable::iterator ci = c1; ci!=c2; ci++) {
        if(keyp((ci->second)->GetKey()))
            ctors.insert(ci->second);
    }
}

void ChannelServer::FindElements (const std::string& id_pat, std::vector<ModelElement*>& elements)
{
    theElements->Find(id_pat,elements);
    sort(elements.begin(),elements.end(),TYPE_CMP());
}

