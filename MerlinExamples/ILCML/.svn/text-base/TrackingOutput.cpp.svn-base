#include "TrackingOutput.h"
#include <fstream>
#include "BasicTransport/NormalTransform.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
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

#define WRITE_FOS(w,p,data) (*fosptr)<<scientific<<setw(w)<<setprecision(p)<<(data)
};


void TrackingOutput::Record(const ComponentFrame* frame, const Bunch* bunch)
{
    if(frame->IsComponent()) {
        string id = (*frame).GetComponent().GetQualifiedName();
		Record(id,bunch);
    }
}

void TrackingOutput::Record(const string& id, const Bunch* bunch)
{
	using std::setw;
	using std::setprecision;

	PSmoments S;
	bunch->GetMoments(S);

	double p0 = bunch->GetReferenceMomentum();
	p0*=1+S.mean(ps_DP);

	double gamma = p0/MeV/ElectronMassMeV;
	double gex = gamma*ProjectedEmittance(S,ps_X,ps_XP);
	double gey = gamma*ProjectedEmittance(S,ps_Y,ps_YP);
	double geyc = S.var(ps_DP) != 0 ? gamma*DispersionCorrectedEmittance(S) : gey;
	double z = bunch->GetReferenceTime();

    double eta_y = S(ps_Y,ps_DP)/S.var(ps_DP);
    double eta_yp = S(ps_YP,ps_DP)/S.var(ps_DP);

	(*fosptr)<<setw(20)<<left<<id<<right;
	WRITE_FOS(17,8,z/1000.0);
	WRITE_FOS(17,8,p0);
	WRITE_FOS(17,8,S.mean(ps_Y)/micrometer);
	WRITE_FOS(17,8,S.std(ps_Y)/micrometer);
	WRITE_FOS(17,8,gey/nanometer);
	WRITE_FOS(17,8,S.std(ps_DP));
	(*fosptr)<<endl;
}

bool TrackingOutput::NewFile(const std::string& fname)
{
	if(fosptr!=0)
		delete fosptr;
	fosptr = new ofstream(fname.c_str());
	return *fosptr;
}



