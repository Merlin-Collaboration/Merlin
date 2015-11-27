#ifndef _h_TrackingOutputASCII
#define _h_TrackingOutputASCII
#include "BeamDynamics/TrackingSimulation.h"
#include <fstream>

// If USE_BOOST_GZ then using a filename ending in ".gz" gives a gzipped file.
#ifdef USE_BOOST_GZ
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif



class TrackingOutputASCII : public SimulationOutput {
	public: 
		TrackingOutputASCII(const std::string& filename);
		~TrackingOutputASCII();
		size_t turn_number;
		void SuppressUnscattered(const int factor);
	
	protected:
		void Record(const ComponentFrame* frame, const Bunch* bunch);
		void RecordInitialBunch(const Bunch* bunch);
		void RecordFinalBunch(const Bunch* bunch);

	private:
#ifdef USE_BOOST_GZ
		ofstream* backing_file;	
		boost::iostreams::filtering_ostream* output_file;
#else
		std::ofstream* output_file;
#endif
		int suppress_factor;
	
};

#endif
