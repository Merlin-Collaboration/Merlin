#include "BeamDynamics/TrackingOutputAV.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

using namespace std;
using namespace ParticleTracking;

TrackingOutputAV::TrackingOutputAV(const std::string& filename):
SimulationOutput(),turn_number(1), suppress_factor(1), current_s(0), current_s_set(0)
{
	output_file = new std::ofstream(filename.c_str());
	(*output_file) << "#id turn S x xp y yp dp type" << std::endl;
}

TrackingOutputAV::~TrackingOutputAV(){
	output_file->close();
	delete output_file;
}

void TrackingOutputAV::Record(const ComponentFrame* frame, const Bunch* bunch){
	if(current_s_set != 1){current_s = frame->GetPosition(); current_s_set = 1; cout << "Setting S = " << current_s << endl; }
	if(frame->GetPosition() < current_s){ turn_number++; cout << "Incrementing turn number = " << turn_number << endl;}
	current_s = frame->GetPosition(); 
	
    //~ if(!frame->IsComponent()) {
    if( !frame->IsComponent() || (turn_number != turn) || (frame->GetPosition() >= end_s || frame->GetPosition() <= start_s) ) {
    //~ if( !frame->IsComponent() || (turn_number != turn) ) {
		return;
	}
	string id = (*frame).GetComponent().GetQualifiedName();
	double zComponent=frame->GetPosition() + frame->GetGeometryLength()/2;

	const ParticleBunch* PB = static_cast<const ParticleBunch*>(bunch);
	for(ParticleBunch::const_iterator pb = PB->begin(); pb!= PB->end(); pb++){
		if(pb->type() != -1){
			(*output_file) 	<< int(pb->id()) << " "
							<< turn_number << " "
							<< std::fixed 
							<< zComponent << " "
							//<< id << " "
							<< std::scientific
							<< pb->x()*1e3  << " "
							<< pb->xp()*1e3  << " "
							<< pb->y()*1e3  << " "
							<< pb->yp()*1e3  << " "
							<< pb->dp()  << " "
							<< std::fixed
							<< int(pb->type()) <<  " "
							<< int(pb->id())  << endl;
		}
	}

}
void TrackingOutputAV::RecordInitialBunch(const Bunch* bunch){}
void TrackingOutputAV::RecordFinalBunch(const Bunch* bunch){}

void TrackingOutputAV::SuppressUnscattered(const int factor){
	suppress_factor = factor;
}
