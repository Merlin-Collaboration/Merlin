#include "DFSOutput.h"
#include "ILCDFS_IO.h"
#include "BasicTransport/NormalTransform.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include <fstream>
#include <iomanip>

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;

namespace {

	double DispersionCorrectedEmittance(const PSmoments& S)
	{
		double s36 = S(ps_Y,ps_DP);
		double s46 = S(ps_YP,ps_DP);
		double dp2 = S.var(ps_DP);
		double s33 = S.var(ps_Y)-s36*s36/dp2;
		double s34 = S(ps_Y,ps_YP)-s36*s46/dp2;
		double s44 = S.var(ps_YP)-s46*s46/dp2;
		return sqrt(s33*s44-s34*s34);
	}

};

#define WRITE_FOS(w,p,data) (*fos)<<scientific<<setw(w)<<setprecision(p)<<(data)

void DFSOutput::RecordInitialBunch(const Bunch* bunch)
{
	Record("INITIAL",bunch);
}

void DFSOutput::RecordFinalBunch(const Bunch* bunch)
{
	Record("FINAL",bunch);
}

void DFSOutput::Record(const ComponentFrame* frm, const Bunch* bunch)
{
	Record((*frm).GetComponent().GetQualifiedName(),bunch);
}

void DFSOutput::Record(const string& id, const Bunch* bunch)
{
	if(!fos)
		return;

	using std::setw;
	using std::setprecision;

	PSmoments S;
	bunch->GetMoments(S);

	double p0 = bunch->GetReferenceMomentum();
	p0*=1+S.mean(ps_DP);

	double gamma = p0/MeV/ElectronMassMeV;
	double gex = gamma*ProjectedEmittance(S,ps_X,ps_XP);
	double gey = gamma*ProjectedEmittance(S,ps_Y,ps_YP);
	double geyc = gamma*DispersionCorrectedEmittance(S);
	double z = bunch->GetReferenceTime();

	(*fos)<<setw(20)<<left<<id<<right;
	WRITE_FOS(12,3,z);
	WRITE_FOS(12,3,p0);
	WRITE_FOS(12,3,gex);
	WRITE_FOS(12,3,gey);
	WRITE_FOS(12,3,geyc);
	WRITE_FOS(12,3,S.mean(ps_Y));
	(*fos)<<endl;
}

void DFSOutput::NewFile(const std::string& fname)
{
	if(fos)
		delete fos;
	fos = new ofstream(fname.c_str());
	if(!*fos) {
		dfs_trace(dfs_trace::error)<<"**** DFSOutput failed to open file "<<fname<<endl;
		delete fos;
		fos=0;
	}
}


