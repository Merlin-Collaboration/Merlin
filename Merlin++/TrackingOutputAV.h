/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_TrackingOutputAV
#define _h_TrackingOutputAV
#include "TrackingSimulation.h"
#include <fstream>

class TrackingOutputAV: public SimulationOutput
{
public:
	TrackingOutputAV(const std::string& filename);
	~TrackingOutputAV();

	/**
	 * This allows us to refrain from outputting non scattered particles
	 */
	void SuppressUnscattered(const bool s);

	/**
	 * This chooses the s coordinate range in which we want to output particle tracks
	 */
	void SetSRange(double start, double end)
	{
		start_s = start;
		end_s = end;
		s_range_set = 1;
	}

	/**
	 * This chooses the turn in which we want to output particle tracks
	 */
	void SetTurn(int inturn)
	{
		turn = inturn;
		single_turn = 1;
	}

	/**
	 * This chooses the turns between which we want to output particle tracks
	 */
	void SetTurnRange(int start, int end)
	{
		start_turn = start;
		end_turn = end;
		turn_range_set = 1;
		single_turn = 0;
	}

protected:
	void Record(const ComponentFrame* frame, const Bunch* bunch);
	void RecordInitialBunch(const Bunch* bunch);
	void RecordFinalBunch(const Bunch* bunch);

private:
	std::ofstream* output_file;

	unsigned int turn;
	double start_s;
	double end_s;
	double current_s;

	bool suppress_unscattered;
	bool current_s_set;
	bool single_turn;
	bool turn_range_set;
	bool s_range_set;

	unsigned int turn_number;
	unsigned int start_turn;
	unsigned int end_turn;

};

#endif
