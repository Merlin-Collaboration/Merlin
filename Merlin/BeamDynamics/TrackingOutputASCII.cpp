#include "TrackingOutputASCII.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

using namespace std;
using namespace ParticleTracking;

TrackingOutputASCII::TrackingOutputASCII(const std::string& filename):
SimulationOutput(),turn_number(0), suppress_factor(1){
	
#ifdef USE_BOOST_GZ				
	backing_file = new std::ofstream(filename.c_str());
	output_file = new boost::iostreams::filtering_ostream;
	if (filename.length()>3 and (0 == filename.compare (filename.length() - 3, 3, ".gz"))){
		output_file->push(boost::iostreams::gzip_compressor(boost::iostreams::zlib::best_compression));
	}
	output_file->push(*backing_file);
#else
	output_file = new std::ofstream(filename.c_str());
#endif
	(*output_file) << "#id turn S x xp y yp dp type" << std::endl;
}


TrackingOutputASCII::~TrackingOutputASCII(){
#ifdef USE_BOOST_GZ
	output_file->strict_sync();
	output_file->pop();
	backing_file->close();
#else
	output_file->close();
#endif
	delete output_file;
}

void TrackingOutputASCII::Record(const ComponentFrame* frame, const Bunch* bunch){
    if(!frame->IsComponent()) {
		return;
	}
	string id = (*frame).GetComponent().GetQualifiedName();
	double zComponent=frame->GetPosition() + frame->GetGeometryLength()/2;

	const ParticleBunch* PB = static_cast<const ParticleBunch*>(bunch);
	for(ParticleBunch::const_iterator pb = PB->begin(); pb!= PB->end(); pb++){
		//~ if (int(pb->type()) < 0 && int(pb->id())%suppress_factor != 0) continue;
		
		(*output_file) << int(pb->id()) << " "
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
							<< int(pb->type())  << endl;
	}

}
void TrackingOutputASCII::RecordInitialBunch(const Bunch* bunch){}
void TrackingOutputASCII::RecordFinalBunch(const Bunch* bunch){}

void TrackingOutputASCII::SuppressUnscattered(const int factor){
	suppress_factor = factor;
}
