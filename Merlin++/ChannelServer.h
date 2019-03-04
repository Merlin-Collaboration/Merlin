/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ChannelServer_h
#define ChannelServer_h 1

#include "merlin_config.h"
#include <string>
#include <map>
#include "AcceleratorModel.h"
#include "ElementRepository.h"

class ModelElement;
class StringPattern;
class RWChannel;
class ROChannel;

/**
 *	Responsible for providing an interface to ROChannel and
 *	RWChannel for the ModelElements in the AcceleratorModel.
 */

class ChannelServer
{
public:
	/**
	 *	Abstract factory class for RWChannel objects.
	 */
	class ChannelCtor
	{
	public:
		/**
		 * destructor
		 */
		virtual ~ChannelCtor()
		{
		}

		/**
		 *	Constructs a channel for the specified ModelElement.
		 */
		virtual ROChannel* ConstructRO(ModelElement* anElement) = 0;

		/**
		 *	Constructs a channel for the specified ModelElement.
		 */
		virtual RWChannel* ConstructRW(ModelElement* anElement) = 0;

		/**
		 *	Returns the ID of the channel (i.e. type.key).
		 *	@return Channel ID
		 */
		std::string GetID();

		/**
		 *	Returns the ModelElement type.
		 *	@return ModelElement type
		 */
		const string& GetType() const;

		/**
		 *	Returns the channel key.
		 *	@return Channel key
		 */
		const string& GetKey() const;

	protected:

		ChannelCtor(const string& aType, const string& aKey);

		std::string type;
		std::string key;
	};

	typedef std::map<std::string, ChannelCtor*> CtorTable;

	/**
	 *	Returns in channels all ROChannels matching chID.
	 *	Returns the number of channels found.
	 *
	 *	@param[out] channels All ROChannels matching chID
	 *	@return Number of channels found
	 */
	size_t GetROChannels(const string& chID, std::vector<ROChannel*>& channels);

	/**
	 *	Returns in channels all RWChannels matching chID.
	 *	Returns the number of channels found.
	 *
	 *	@param[out] channels All RWChannels matching chID
	 *	@return Number of channels found
	 */
	size_t GetRWChannels(const string& chID, std::vector<RWChannel*>& channels);

	/**
	 *	Returns read-only channels matching chid for all
	 *	matching components in aBeamline. Note that only
	 *	channels associated with AcceleratorComponents can be
	 *	extracted using this method.
	 *
	 *	@param[out] channels All ROChannels matching chID for all matching
	 *	components in aBeamline
	 */
	size_t GetROChannels(AcceleratorModel::Beamline& aBeamline, const std::string& chid,
		std::vector<ROChannel*>& channels);

	/**
	 *	Returns read-write channels matching chid for all
	 *	matching components in aBeamline. Note that only
	 *	channels associated with AcceleratorComponents can be
	 *	extracted using this method.
	 *	@param[out] channels All ROChannels matching chID for all matching
	 *	components in aBeamline
	 */
	size_t GetRWChannels(AcceleratorModel::Beamline& aBeamline, const std::string& chid,
		std::vector<RWChannel*>& channels);

	/**
	 *	Adds a ChannelCtor object to the server.
	 */
	void RegisterCtor(ChannelCtor* chctor);

	void SetRepository(ElementRepository* me_repo);

	~ChannelServer();

private:

	ElementRepository* theElements;
	CtorTable chCtors;

	void FindCtors(const string& type, const string& keypat, std::set<ChannelCtor*>& ctors);

	/**
	 *	Returns in elements the ModelElements that match
	 *	pattern. elements is sorted by TYPE.
	 *
	 *	@param[out] elements ModelElements that match pattern id_pat, sorted by
	 *	type
	 */
	void FindElements(const std::string& id_pat, std::vector<ModelElement*>& elements);
};

inline ChannelServer::ChannelCtor::ChannelCtor(const string& aType, const string& aKey) :
	type(aType), key(aKey)
{
}

inline const string& ChannelServer::ChannelCtor::GetType() const
{
	return type;
}

inline const string& ChannelServer::ChannelCtor::GetKey() const
{
	return key;
}

#endif
