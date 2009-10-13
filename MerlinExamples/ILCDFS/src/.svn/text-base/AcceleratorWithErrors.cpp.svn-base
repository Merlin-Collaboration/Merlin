/////////////////////////////////////////////////////////////////////////
// Class AcceleratorWithErrors implementation
// Represents the physical accelerator with errors.
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2008/01/14 21:08:22 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorWithErrors.h"
#include "utility/StringPattern.h"
#include "AcceleratorModel/Frames/LatticeFrame.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "ILCDFS_IO.h"
#include "Random/RandomNG.h"
#include "AcceleratorModel/ActiveMonitors/BPM.h"

class AlignmentError : public FrameTraverser {
public:
	AlignmentError(const string& sp, bool clear_t)
		: fpattern(sp),s(sp),ct(clear_t) {}
	void ActOn(LatticeFrame* f);

	void Init() { count=0; }
	void Report() const;

protected:
	virtual void ApplyError(LatticeFrame* f)=0;
	size_t count;
private:
	StringPattern fpattern;
	string s;
	bool ct;
};

// implementation
void AlignmentError::ActOn(LatticeFrame* f)
{
	const string& id = f->GetQualifiedName();
	if(fpattern(id)) {
		dfs_trace(dfs_trace::level_3)<<"adding alignment errors to "<<id<<endl;
		if(ct)
			f->ClearLocalFrameTransform();
		ApplyError(f);
		count++;
	}
}

void AlignmentError::Report() const
{
	dfs_trace(dfs_trace::level_2)<<count<<" alignment errors added to "<<s<<endl;
}

namespace {

	class TransverseError : public AlignmentError {
	public:
		TransverseError(const string& sp, double xrms, double yrms, bool ct)
			: AlignmentError(sp,ct),xvar(xrms*xrms),yvar(yrms*yrms) {}
	protected:
		void ApplyError(LatticeFrame* f);
	private:
		double xvar, yvar;
	};

	class RotationError : public AlignmentError {
	public:
		RotationError(const string& sp, double xrms, double yrms, double zrms, bool ct)
			: AlignmentError(sp,ct),xvar(xrms*xrms),yvar(yrms*yrms),zvar(zrms*zrms) {}
	protected:
		void ApplyError(LatticeFrame* f);
	private:
		double xvar, yvar, zvar;
	};

	void TransverseError::ApplyError(LatticeFrame* f)
	{
		double x = RandomNG::normal(0,xvar);
		double y = RandomNG::normal(0,yvar);
		f->Translate(x,y,0);
	}

	void RotationError::ApplyError(LatticeFrame* f)
	{
		double x = RandomNG::normal(0,xvar);
		double y = RandomNG::normal(0,yvar);
		double z = RandomNG::normal(0,zvar);
		f->RotateX(x);
		f->RotateY(y);
		f->RotateZ(z);
	}
};

AcceleratorWithErrors::AcceleratorWithErrors(const std::string& name, AcceleratorModel* mdl, BeamData* beam0)
:Accelerator(name,mdl,beam0)
{}

void AcceleratorWithErrors::TransverseErrors(const string& pattern, double x_rms, double y_rms, bool clear_transform)
{
	alignErrorDefs.push_back(new TransverseError(pattern,x_rms,y_rms,clear_transform));
}

void AcceleratorWithErrors::RotationErrors(const string& pattern, double x_rms, double y_rms, double z_rms, bool clear_transform)
{
	alignErrorDefs.push_back(new RotationError(pattern,x_rms,y_rms,z_rms,clear_transform));
}

void AcceleratorWithErrors::MagnetScaleError(const string& pattern, double rms)
{
	dfs_trace(dfs_trace::warning)<<"WARNING! *** Magnet scale errors are not yet implemented"<<endl;
}

void AcceleratorWithErrors::KlystronErrors(const string& pattern, double v_rms, double phi_rms)
{
	dfs_trace(dfs_trace::warning)<<"WARNING! *** Klystron errors are not yet implemented"<<endl;
}

void AcceleratorWithErrors::BPMresolution(double rms)
{
	dfs_trace(dfs_trace::level_1)<<"BPM resolution = "<<rms<<endl;
	AcceleratorModel::Beamline bl = itsAccModel->GetBeamline();
	for(AcceleratorModel::BeamlineIterator c = bl.begin(); c!=bl.end(); c++) {
		if((*c)->IsComponent()) {
			BPM* bpm = dynamic_cast<BPM*>(&((*c)->GetComponent()));
			if(bpm)
				bpm->SetRes(rms);
		}
	}
}

void AcceleratorWithErrors::BPMlinearError(double rms)
{
	bpmCalErr2 = rms*rms;
}

void AcceleratorWithErrors::InitialBeamJitterInSigma(double xrms, double yrms)
{
	// TO DO
	dfs_trace(dfs_trace::warning)<<"WARNING! *** beam jitter is not yet implemented"<<endl;
}

void AcceleratorWithErrors::ApplyStaticErrors()
{
	for(list<AlignmentError*>::iterator ei = alignErrorDefs.begin(); ei!=alignErrorDefs.end(); ei++) {
		AlignmentError* aerr = *ei;
		aerr->Init();
		(itsAccModel->GetGlobalFrame()).Traverse(*aerr);
		aerr->Report();
	}

	// Add BPM calibration errors, if applicable
	dfs_trace(dfs_trace::level_2)<<"BPM rms calibration error = "<<sqrt(bpmCalErr2)<<endl;


	AcceleratorModel::Beamline bl = itsAccModel->GetBeamline();
	for(AcceleratorModel::BeamlineIterator c = bl.begin(); c!=bl.end(); c++) {
		if((*c)->IsComponent()) {
			BPM* bpm = dynamic_cast<BPM*>(&((*c)->GetComponent()));
			if(bpm) {
				double ce = bpmCalErr2==0? 1.0 : 1.0+RandomNG::normal(0,bpmCalErr2);
				bpm->SetScale(ce,ce); // Same error in both planes?
			}
		}
	}
}
