#include <iostream>
#include <sstream>
#include <string>

#include "AcceleratorModel/ApertureSurvey.h"
#include "AcceleratorModel/StdComponent/Collimator.h"

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

using namespace std;

ApertureSurvey::ApertureSurvey(AcceleratorModel* model, string directory, double step, size_t points_per_element)
{
	step_size = step;
	points = points_per_element;
	AccMod = model;
	
	//Create a file in the given directory
	ostringstream file_stream_name;
	file_stream_name << directory <<"Aperture_Survey_"<< step << "_steps_OR_" << points << "_points.txt";	
	ofstream* output_file = new ofstream(file_stream_name.str().c_str());	
	if(!output_file->good())    {
        std::cerr << "Could not open ApertureSurvey file" << std::endl;
        exit(EXIT_FAILURE);
    }   
	(*output_file) << "#name\ttype\ts_end\tlength\tap_px\tap_mx\tap_py\tap_my" << endl;
	
	//~ DoApertureSurvey();
	double s = 0;
	double last_sample = 0-step_size;
	double lims[4];
	//~ cout << "aperture_survey" << endl;
		
	for (AcceleratorModel::BeamlineIterator bi = AccMod->GetBeamline().begin(); bi != AccMod->GetBeamline().end(); bi++){
		
		AcceleratorComponent *ac = &(*bi)->GetComponent();

		if (fabs(s - ac->GetComponentLatticePosition())> 1e-6){exit(1);}

		Aperture* ap =	ac->GetAperture();
		Collimator* aCollimator = dynamic_cast<Collimator*>(ac);
		
		//cout << "ap "<< s << " " << ac->GetName() << " " << ac->GetLength() << ((aCollimator!=NULL)?" Collimator":"");
		//if (ap != NULL) cout << " " << ap->GetMaterial();
		//cout << endl;
		
		vector<double> zs;
		if (points > 0){
			for (size_t i=0; i<points; i++) {zs.push_back(i*ac->GetLength()*(1.0/(points-1)));}
		}
		else{
			//~ cout << "last_sample = " << last_sample <<endl;
			while (last_sample+step_size < s + ac->GetLength()){
				last_sample += step_size;
				//~ //cout << "add step " <<  last_sample << endl;
				zs.push_back(last_sample-s);
			}	
		}
		for (size_t zi = 0; zi < zs.size(); zi++){
			double z = zs[zi];
			//~ cout << "call check_aperture(" << z+s << ")" << endl;
			if (ap != NULL) {
				ApertureSurvey::CheckAperture(ap, z, lims);
			}
			else{
				lims[0] = lims[1] = lims[2]= lims[3] = 1;
			}
			(*output_file) << ac->GetName() << "\t";
			(*output_file) << ac->GetType() << "\t";
			(*output_file) << ac->GetComponentLatticePosition()+ac->GetLength() << "\t";
			(*output_file) << ac->GetLength() << "\t";
			//~ (*os) << s+z << "\t";
			(*output_file) << lims[0] << "\t";
			(*output_file) << lims[1] << "\t";
			(*output_file) << lims[2] << "\t";
			(*output_file) << lims[3] << endl;
		}
		s += ac->GetLength();
	}
}

ApertureSurvey::ApertureSurvey(AcceleratorModel* model, std::ostream* os, double step, size_t points_per_element)
{
	step_size = step;
	points = points_per_element;
	AccMod = model;
	
	(*os) << "#name\ttype\ts_end\tlength\tap_px\tap_mx\tap_py\tap_my" << endl;
	
	//~ DoApertureSurvey();
	double s = 0;
	double last_sample = 0-step_size;
	double lims[4];
	//~ cout << "aperture_survey" << endl;
		
	for (AcceleratorModel::BeamlineIterator bi = AccMod->GetBeamline().begin(); bi != AccMod->GetBeamline().end(); bi++){
		
		AcceleratorComponent *ac = &(*bi)->GetComponent();

		if (fabs(s - ac->GetComponentLatticePosition())> 1e-6){exit(1);}

		Aperture* ap =	ac->GetAperture();
		Collimator* aCollimator = dynamic_cast<Collimator*>(ac);
		
		//cout << "ap "<< s << " " << ac->GetName() << " " << ac->GetLength() << ((aCollimator!=NULL)?" Collimator":"");
		//if (ap != NULL) cout << " " << ap->GetMaterial();
		//cout << endl;
		
		vector<double> zs;
		if (points > 0){
			for (size_t i=0; i<points; i++) {zs.push_back(i*ac->GetLength()*(1.0/(points-1)));}
		}
		else{
			//~ cout << "last_sample = " << last_sample <<endl;
			while (last_sample+step_size < s + ac->GetLength()){
				last_sample += step_size;
				//~ //cout << "add step " <<  last_sample << endl;
				zs.push_back(last_sample-s);
			}	
		}
		for (size_t zi = 0; zi < zs.size(); zi++){
			double z = zs[zi];
			//~ cout << "call check_aperture(" << z+s << ")" << endl;
			if (ap != NULL) {
				ApertureSurvey::CheckAperture(ap, z, lims);
			}
			else{
				lims[0] = lims[1] = lims[2]= lims[3] = 1;
			}
			
			(*os) << ac->GetName() << "\t";
			(*os) << ac->GetType() << "\t";
			(*os) << ac->GetComponentLatticePosition()+ac->GetLength() << "\t";
			(*os) << ac->GetLength() << "\t";
			//~ (*os) << s+z << "\t";
			(*os) << lims[0] << "\t";
			(*os) << lims[1] << "\t";
			(*os) << lims[2] << "\t";
			(*os) << lims[3] << endl;
		}
		s += ac->GetLength();
	}
}

void ApertureSurvey::CheckAperture(Aperture* ap, double s, double *aps){
	//~ cout << "CheckAperture" << endl;
	const double step = 1e-6;
	const double max = 1.0;
	const double min = 0.0;
	
	// iterate through directions
	for (int dir=0; dir<4; dir++){
		double xdir=0, ydir=0;
		if (dir==0) {xdir=+1;}
		else if (dir==1) {xdir=-1;}
		else if (dir==2) {ydir=+1;}
		else if (dir==3) {ydir=-1;}
		
		aps[dir] = 0;
		
		// scan for limit
		double below=min, above=max;
		
		while(above-below > step){
			double guess = (above+below)/2;
			
			if (ap->PointInside(xdir*guess, ydir*guess, s)){
				below = guess;
			}
			else{
				above = guess;
			}
		}
		aps[dir] = (above+below)/2;
	}
}
