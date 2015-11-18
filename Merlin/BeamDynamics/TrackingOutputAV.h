#ifndef _h_TrackingOutputAV
#define _h_TrackingOutputAV
#include "BeamDynamics/TrackingSimulation.h"
#include <fstream>

class TrackingOutputAV : public SimulationOutput {
	public: 
		TrackingOutputAV(const std::string& filename);
		~TrackingOutputAV();
		void SuppressUnscattered(const int factor);
		
		// This chooses the s coordinate range in which we want to output particle tracks
		void SetSRange(double start, double end){start_s = start; end_s = end;}
		
		// This chooses the turn in which we want to output particle tracks
		void SetTurn(int inturn){turn = inturn;}	
	
	protected:
		void Record(const ComponentFrame* frame, const Bunch* bunch);
		void RecordInitialBunch(const Bunch* bunch);
		void RecordFinalBunch(const Bunch* bunch);

	private:
		std::ofstream* output_file;

		int suppress_factor;
		int turn;
		double start_s;
		double end_s;
		double current_s;
		bool current_s_set;
		int turn_number;
	
};

#endif
