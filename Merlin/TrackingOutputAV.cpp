/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "TrackingOutputAV.h"
#include "ParticleBunch.h"

using namespace std;
using namespace ParticleTracking;

TrackingOutputAV::TrackingOutputAV(const std::string& filename) :
	SimulationOutput(), current_s(0), suppress_unscattered(1), current_s_set(0), single_turn(1), turn_range_set(0),
	s_range_set(0), turn_number(0)
{
	output_file = new std::ofstream(filename.c_str());
	(*output_file) << "#id turn S x xp y yp dp type" << std::endl;

	// This sets a flag used in the TrackingSimulation class
	output_all = 1;
}

TrackingOutputAV::~TrackingOutputAV()
{
	output_file->close();
	delete output_file;
}

void TrackingOutputAV::Record(const ComponentFrame* frame, const Bunch* bunch)
{

	// if no starting s value is set we store it
	if(current_s_set != 1)
	{
		current_s = frame->GetPosition();
		turn_number = 1;
		current_s_set = 1;
		cout << "\n\nSetting S = " << current_s << endl;
	}

	// if the position of the element is less than the stored s position we increment the turn, and continue to store s
	if(frame->GetPosition() < current_s)
	{
		turn_number++;
		cout << "Incrementing turn number = " << turn_number << endl;
	}
	current_s = frame->GetPosition();

	// if the current frame is not a component we exit
	if(!frame->IsComponent())
	{
		return;
	}

	// check if we are using a single turn & on the specified turn (default is single turn, with turn = 1)
	//~ if( single_turn && (turn_number != turn) ){return;}

	// check if we are using a turn range & we are within the range
	if(turn_range_set && (turn_number > end_turn || turn_number < start_turn))
	{
		return;
	}

	// check if we are using an s range & we are within the range
	if(s_range_set && (frame->GetPosition() < start_s || frame->GetPosition() > end_s))
	{
		return;
	}

	string id = (*frame).GetComponent().GetQualifiedName();
	double zComponent = frame->GetPosition() + frame->GetGeometryLength() / 2;

	const ParticleBunch* PB = static_cast<const ParticleBunch*>(bunch);
	for(ParticleBunch::const_iterator pb = PB->begin(); pb != PB->end(); pb++)
	{
		if((suppress_unscattered && pb->type() != -1) || (!suppress_unscattered))
		{
			(*output_file) << int(pb->id()) << " "
						   << turn_number << " "
						   << std::fixed
						   << zComponent << " "
			    //<< id << " "
						   << std::scientific
						   << pb->x() * 1e3  << " "
						   << pb->xp() * 1e3  << " "
						   << pb->y() * 1e3  << " "
						   << pb->yp() * 1e3  << " "
						   << pb->dp()  << " "
						   << std::fixed
						   << int(pb->type()) <<  " "
						   << int(pb->id())  << endl;
		}
	}

}

void TrackingOutputAV::RecordInitialBunch(const Bunch* bunch)
{
}
void TrackingOutputAV::RecordFinalBunch(const Bunch* bunch)
{
}

void TrackingOutputAV::SuppressUnscattered(const bool s)
{
	suppress_unscattered = s;
}
