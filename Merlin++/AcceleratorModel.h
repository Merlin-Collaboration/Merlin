/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef AcceleratorModel_h
#define AcceleratorModel_h 1

#include "merlin_config.h"
#include <vector>
#include <string>
#include "ComponentFrame.h"

#include "AcceleratorComponent.h"
#include "ModelElement.h"
#include "ElementRepository.h"
#include "LatticeFrame.h"
#include "StringPattern.h"
#include "ring_iterator.h"
#include "MerlinException.h"
#include "AcceleratorSupport.h"

class ChannelServer;
class ComponentFrame;
template<class T> class TComponentFrame;
class RWChannel;
class ROChannel;

class AcceleratorModel
{
public:

	typedef size_t Index;

	/**
	 * A sequence of ComponentFrame objects representing the
	 * complete accelerator lattice.
	 */
	typedef std::vector<ComponentFrame*> FlatLattice;

	// Iterator definitions.
	typedef FlatLattice::iterator BeamlineIterator;
	typedef FlatLattice::const_iterator ConstBeamlineIterator;

	/**
	 * Represents the complete or contiguous sub-section of
	 * the accelerator lattice.
	 */
	class Beamline
	{
	protected:

		/**
		 * Wrapper used by for template<T> void Track(T&) function.
		 */
		template<class T>
		class TRK
		{
		public:
			explicit TRK(T& aT) :
				_t(aT)
			{
			}

			void operator ()(ComponentFrame* frame)
			{
				_t(frame);
			}
			T& _t;
		};

	public:
		Beamline(BeamlineIterator fst, BeamlineIterator lst, Index fst_i, Index lst_i) :
			first(fst), last(lst), first_i(fst_i), last_i(lst_i)
		{
		}

		Beamline() :
			first(), last(), first_i(), last_i()
		{
		}

		/**
		 * Checks if the beamline is reversed.
		 * @return Returns a bool set to true if the beamline is reversed.
		 */
		bool IsReversed()
		{
			return first > last;
		}

		/**
		 * Template function that iterates the functor object tobj
		 * over the Beamline. Returns a reference to tobj on exit.
		 */
		template<class T> T& Track(T& tobj)
		{
			/**
			 * We make use of the TRK wrapper, to avoid the call to tobj's
			 * copy constructor
			 */
			std::for_each(begin(), end(), TRK<T>(tobj));
			return tobj;
		}

		/**
		 * Standard library type iterator accessors.
		 */

		/**
		 * Gets a BeamlineIterator pointing to the first element in the Beamline.
		 * @return A BeamlineIterator pointing to the first element in the Beamline.
		 */
		BeamlineIterator begin()
		{
			return first;
		}

		/**
		 * Gets a const BeamlineIterator pointing to the first element in the Beamline.
		 * @return A const BeamlineIterator pointing to the first element in the Beamline.
		 */
		ConstBeamlineIterator begin() const
		{
			return first;
		}

		/**
		 * Gets a BeamlineIterator pointing to the end of the Beamline.
		 * @return A BeamlineIterator pointing to the end of the Beamline.
		 */
		BeamlineIterator end()
		{
			return last + 1;
		}

		/**
		 * Gets a const BeamlineIterator pointing to the end of the Beamline.
		 * @return A const BeamlineIterator pointing to the end of the Beamline.
		 */
		ConstBeamlineIterator end() const
		{
			return last + 1;
		}

		/**
		 * Gets a reference to the first ComponentFrame.
		 * @return A reference to the first ComponentFrame.
		 */
		ComponentFrame& First()
		{
			return **first;
		}

		/**
		 * Gets a const reference to the first ComponentFrame.
		 * @return A const reference to the first ComponentFrame.
		 */
		const ComponentFrame& First() const
		{
			return **first;
		}

		/**
		 * Gets a reference to the last ComponentFrame.
		 * @return A reference to the last ComponentFrame.
		 */
		ComponentFrame& Last()
		{
			return **last;
		}

		/**
		 * Gets a const reference to the last ComponentFrame.
		 * @return A const reference to the last ComponentFrame.
		 */
		const ComponentFrame& Last() const
		{
			return **last;
		}

		/**
		 * Gets the index of the first element in the Beamline.
		 * @return The index of the first element in the Beamline.
		 */
		Index first_index() const
		{
			return first_i;
		}

		/**
		 * Gets the index of the last element in the Beamline.
		 * @return The index of the last element in the Beamline.
		 */
		Index last_index() const
		{
			return last_i;
		}

	private:
		// Data Members for Class Attributes
		BeamlineIterator first;

		BeamlineIterator last;

		Index first_i;
		Index last_i;

	}; //End Beamline class

	/**
	 * Exception class for GetBeamline functions
	 */
	class BadRange: public MerlinException
	{
	public:
		BadRange() :
			MerlinException("Bad Beamline Range")
		{
		}
	};

	typedef ring_iterator<FlatLattice, FlatLattice::iterator> RingIterator;

	/**
	 * Constructor.
	 */
	AcceleratorModel();

	/**
	 * Destructor.
	 */
	~AcceleratorModel();

	/**
	 * Returns the entire beamline of the model.
	 * @return The entire Beamline for this AcceleratorModel.
	 */
	Beamline GetBeamline();

	/**
	 * Returns the beamline from elements n1 to n2.
	 * @param[in] n1 The index of the first element to return in the Beamline.
	 * @param[in] n2 The index of the last element to return in the Beamline.
	 * @exception Throws a BadRange exception if the requested range cannot be found.
	 * @return The requested Beamline.
	 */
	Beamline GetBeamline(Index n1, Index n2);

	/**
	 * Returns a Beamline from the n1-th occurrence of the
	 * component whose qualified name  matches the pattern
	 * pat1, to the n2-th occurrence of the component matching
	 * patl2. Throws BadRange if no section is found.
	 * @param[in] pat1 A string containing the qualified name of the elements to match for the start of the Beamline.
	 * @param[in] pat2 A string containing the qualified name of the elements to match for the end of the Beamline.
	 * @param[in] n1 The number of the occurrence of pat1 from which to start the Beamline (default is the first).
	 * @param[in] n2 The number of the occurrence of pat2 from which to end the Beamline (default is the first).
	 * @exception Throws a BadRange exception if the requested range cannot be found.
	 * @return The requested Beamline.
	 */
	Beamline GetBeamline(const std::string& pat1, const std::string& pat2, int n1 = 1, int n2 = 1);

	/**
	 * Assumes that the AcceleratorModel represents a ring
	 * accelerator, and returns a RingIterator which iterates
	 * continuously the ring. n represents the offset from the
	 * beginning of the ring as defined in the AcceleratorModel.
	 * @param[in] n An integer containing the offset from the beginning of the ring in number of elements.
	 * @return A RingIterator for the entire ring starting at position n.
	 */
	RingIterator GetRing(int n = 0);

	/**
	 * Returns the reversed complete beamline of the model.
	 * @return The reversed Beamline
	 */
	Beamline GetReversedBeamline();

	/**
	 * Returns in results all ComponentFrame objects whose name
	 * matches the string pattern pat. Returns the length of
	 * results on exit. Note that the previous contents of
	 * results is overwritten. Components are returned in
	 * Beamline order.
	 * @param[in] pat The pattern of element names that should be matched.
	 * @param[out] results The vector of ComponentFrame pointers that will contain the output results.
	 * @return An integer containing the number of extracted components.
	 */
	int ExtractComponents(const std::string& pat, std::vector<ComponentFrame*>& results);

	/**
	 * Returns in results all ModelElement objects whose name
	 * matches the string pattern pat. Returns the length of
	 * results on exit. Note that the previous contents of
	 * results is overwritten. The order results is undefined.
	 * @param[in] pat The pattern of element names that should be matched.
	 * @param[out] results The vector of ModelElement pointers that will contain the output results.
	 * @return An integer containing the number of extracted model elements.
	 */
	int ExtractModelElements(const std::string& pat, std::vector<ModelElement*>& results);

	/**
	 * template function returning TComponentFrame objects
	 * corresponding to AcceleratorComponents of type T.
	 * pattern is optional string pattern which can be used to
	 * match only those components with a specific
	 * (unqualified) name. Components are returned in Beamline
	 * order.
	 * @param[out] results A vector container for pointers to the type of AcceleratorComponents one wishes to match.
	 * @param[in] pattern A string containing the pattern of unqualified element names to match. Default is to match all elements.
	 * @return An integer containing the number of matched components
	 */
	template<class T> int ExtractTypedComponents(std::vector<TComponentFrame<T>*>& results, const std::string& pattern =
		"*")
	{
		StringPattern p(pattern);
		for(BeamlineIterator i = lattice.begin(); i != lattice.end(); i++)
		{
			TComponentFrame<T>* cf = dynamic_cast<TComponentFrame<T>*>(*i);
			if(cf && p(cf->GetName()))
			{
				results.push_back(cf);
			}
		}
		return results.size();
	}

	/**
	 * template function returning ModelElements of type T.
	 * pattern is optional string pattern which can be used to
	 * match only those components with a specific
	 * (unqualified) name. Order is undefined.
	 * @param[out] results A vector container for pointers to the type of ModelElements one wishes to match.
	 * @param[in] pattern A string containing the pattern of unqualified element names to match. Default is to match all elements.
	 * @return An integer containing the number of matched components
	 */
	template<class T> int ExtractTypedElements(T& results, const std::string& pattern = "*")
	{
		typedef typename T::value_type value_type;
		StringPattern p(pattern);
		for(ElementRepository::iterator i = theElements->begin(); i != theElements->end(); i++)
		{
			value_type mi = dynamic_cast<value_type>(*i);
			if(mi && p(mi->GetName()))
			{
				results.push_back(mi);
			}
		}
		return results.size();
	}

	/**
	 * Returns the indices of components matching pat in iarray
	 * for the entire beamline. iarray is overwritten by this function.
	 * Function returns length of iarray.
	 * @param[in] pat A string containing the pattern of element names to match.
	 * @param[out] iarray A vector array holding the Indices of the element locations that matched the name pattern.
	 * @return A size_t containing the number of matched elements
	 */
	size_t GetIndexes(const std::string& pat, std::vector<Index>& iarray) const;

	/**
	 * Returns the indexes of components matching pat for the entire beamline.
	 * @param[in] pat A string containing the pattern of element names to match.
	 * @return A vector holding the Indices of the element locations that matched the name pattern.
	 */
	std::vector<Index> GetIndexes(const std::string& pat) const;

	/**
	 * Same as above, but limits search to the specified (sub-)beamline.
	 * @param[in] aBeamline The (sub) Beamline to search.
	 * @param[in] pat A string containing the pattern of element names to match.
	 * @param[out] iarray A vector array holding the Indices of the element locations that matched the name pattern.
	 * @return A size_t containing the number of matched elements.
	 */
	size_t GetIndexes(const Beamline& aBeamline, const std::string& pat, std::vector<Index>& iarray) const;

	/**
	 * Returns in channels all ROChannels matching chID.
	 * Returns the number of channels found.
	 * @param[in] chID A string containing the pattern of channel names to match.
	 * @param[out] channels A vector array holding the ROChannels that matched the chID pattern.
	 * @return A size_t containing the number of matched channels.
	 */
	size_t GetROChannels(const std::string& chID, std::vector<ROChannel*>& channels);

	/**
	 * Returns in channels all RWChannels matching chID.
	 * Returns the number of channels found.
	 * @param[in] chID A string containing the pattern of channel names to match.
	 * @param[out] channels A vector array holding the RWChannels that matched the chID pattern.
	 * @return A size_t containing the number of matched channels.
	 */
	size_t GetRWChannels(const std::string& chID, std::vector<RWChannel*>& channels);

	/**
	 * Returns read-only channels matching chid for all
	 * matching components in aBeamline. Note that only
	 * channels associated with AcceleratorComponents can be
	 * extracted using this method.
	 * @param[in] aBeamline The (sub) Beamline to search.
	 * @param[in] chID A string containing the pattern of channel names to match.
	 * @param[out] channels A vector array holding the RWChannels that matched the chID pattern.
	 * @return A size_t containing the number of matched channels.
	 */
	size_t GetROChannels(Beamline& aBeamline, const std::string& chID, std::vector<ROChannel*>& channels);

	/**
	 * Returns read-write channels matching chid for all
	 * matching components in aBeamline. Note that only
	 * channels associated with AcceleratorComponents can be
	 * extracted using this method.
	 * @param[in] aBeamline The (sub) Beamline to search.
	 * @param[in] chID A string containing the pattern of channel names to match.
	 * @param[out] channels A vector array holding the RWChannels that matched the chID pattern.
	 * @return A size_t containing the number of matched channels.
	 */
	size_t GetRWChannels(Beamline& aBeamline, const std::string& chID, std::vector<RWChannel*>& channels);

	/**
	 * Returns the top-level LatticeFrame (global frame) for
	 * the model. The global frame is the root object of the
	 * lattice frame hierarchy.
	 * @return The LatticeFrame containing the globalFrame.
	 */
	LatticeFrame& GetGlobalFrame()
	{
		return *globalFrame;
	}

	/**
	 * Returns the top-level LatticeFrame (global frame) for
	 * the model. The global frame is the root object of the
	 * lattice frame hierarchy.
	 * const version.
	 * @return The LatticeFrame containing the globalFrame.
	 */
	const LatticeFrame& GetGlobalFrame() const
	{
		return *globalFrame;
	}

	/**
	 * Allows clients to construct and add new ModelElement
	 * objects to the AcceleratorModel.
	 * @param[in] element A pointer to the ModelElement to add to the AcceleratorModel.
	 */
	void AddModelElement(ModelElement* element);

	/**
	 * Allows clients to install a new ModelElement into the AcceleratorModel.
	 * The new element position must fall within a Drift, which will be split to make
	 * space.
	 * @param[in] element A pointer to the ModelElement to add to the AcceleratorModel.
	 * @param[in] at the location to place the new element.
	 */
	void InstallModelElement(AcceleratorComponent* element, double at);

	/**
	 * Prints to the specified stream statistics about the
	 * model.
	 * @param[out] os The stream to output the model statistics to.
	 */
	void ReportModelStatistics(std::ostream& os) const;

	/**
	 * Access to AcceleratorSupport objects.
	 * @param[out] supports A reference to the requested output location for the accelerator support list.
	 * @return The number of found accelerator supports found.
	 */
	size_t GetAcceleratorSupports(AcceleratorSupportList& supports);

	/**
	 * Find the lattice position of a given element.
	 * @param[in] RequestedElement The name of the requested element to find.
	 * @return An integer containing the number of the element in the lattice.
	 */
	int FindElementLatticePosition(std::string RequestedElement);

private:

	FlatLattice lattice;
	LatticeFrame* globalFrame;
	ElementRepository* theElements;
	ChannelServer* chServer;

	friend class AcceleratorModelConstructor;

	//Disable copying
	AcceleratorModel(const AcceleratorModel& am);
	AcceleratorModel& operator=(const AcceleratorModel& am);

};

#endif
