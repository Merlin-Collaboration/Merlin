/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/08/17 08:21:59 $
// $Revision: 1.7 $
// 
/////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <vector>
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/RingDeltaTProcess.h"
#include "AcceleratorModel/AcceleratorModel.h"
#include "RingDynamics/ClosedOrbit.h"
#include "RingDynamics/TransferMatrix.h"
#include "NumericalUtils/MatrixPrinter.h"
#include "NumericalUtils/NumericalConstants.h"
#include "TLAS/TLAS.h"
#include "LatticeFunctions.h"

using namespace ParticleTracking;
using namespace TLAS;

LatticeFunction::LatticeFunction(int _i, int _j, int _k) : i(_i), j(_j), k(_k)
{ };

LatticeFunction::~LatticeFunction()
{
    value.clear();
};

void LatticeFunction::GetIndices(int& _i, int& _j, int& _k)
{
    _i = i;
    _j = j;
    _k = k;
};

void LatticeFunction::AppendValue(double v)
{
    value.push_back(v);
};

void LatticeFunction::ClearValues()
{
    value.clear();
};

double LatticeFunction::GetValue(int n)
{
    return value[n];
};

vector<double>::iterator LatticeFunction::begin()
{
    return value.begin();
}

vector<double>::iterator LatticeFunction::end()
{
    return value.end();
}

int LatticeFunction::size()
{
    return value.size();
}

void LatticeFunction::Derivative(LatticeFunction* lfnM, LatticeFunction* lfnP, double dp)
{
    int i, j, k;
    lfnM->GetIndices(i, j, k);
    bool isPathLength = (i==0 && j==0 && k==0) ? true : false;

    vector<double>::iterator lfnMit = lfnM->begin();
    vector<double>::iterator lfnPit = lfnP->begin();

    double dv;
    for(;lfnPit!=lfnP->end();lfnPit++) {

        if(isPathLength)
            dv = *lfnMit;
        else
            dv = (*lfnPit - *lfnMit)/dp;

        value.push_back( dv );
        lfnMit++;
    };
};

LatticeFunctionTable::LatticeFunctionTable(AcceleratorModel* aModel, double refMomentum)
        : theModel(aModel), p0(refMomentum), delta(1.0e-8), bendscale(1.0e-16), symplectify(false), orbitonly(true)
{
    UseDefaultFunctions();
};

struct DeleteLatticeFunction
{
public:
    DeleteLatticeFunction() {};
    void operator()(LatticeFunction* lfn) { delete lfn; };
};

LatticeFunctionTable::~LatticeFunctionTable()
{
    for_each(lfnlist.begin(),lfnlist.end(),DeleteLatticeFunction());
    lfnlist.clear();
};

void LatticeFunctionTable::SetDelta(double new_delta)
{
    delta = new_delta;
};

void LatticeFunctionTable::ScaleBendPathLength(double scale)
{
    bendscale = scale;
};

void LatticeFunctionTable::MakeTMSymplectic(bool flag)
{
    symplectify = flag;
};

void LatticeFunctionTable::AddFunction(int i, int j, int k)
{
    LatticeFunction* lfn = new LatticeFunction(i, j, k);
    lfnlist.push_back(lfn);
    if(k!=0)
        orbitonly = false;
};

void LatticeFunctionTable::UseDefaultFunctions()
{
    RemoveAllFunctions();
    AddFunction(1,0,0); // closed orbit: x
    AddFunction(2,0,0); // closed orbit: px
    AddFunction(3,0,0); // closed orbit: y
    AddFunction(4,0,0); // closed orbit: py
    AddFunction(5,0,0); // closed orbit: ct
    AddFunction(6,0,0); // closed orbit: dp
    AddFunction(1,1,1); // beta_x
    AddFunction(1,2,1); // -alfa_x
    AddFunction(3,3,2); // beta_y
    AddFunction(3,4,2); // -alfa_y
};

void LatticeFunctionTable::UseOrbitFunctions()
{
    RemoveAllFunctions();
    AddFunction(1,0,0); // closed orbit: x
    AddFunction(2,0,0); // closed orbit: px
    AddFunction(3,0,0); // closed orbit: y
    AddFunction(4,0,0); // closed orbit: py
};

void LatticeFunctionTable::Size(int& rows, int& cols)
{
    rows = NumberOfRows();
    cols = lfnlist.size();
};

void LatticeFunctionTable::RemoveAllFunctions()
{
    for_each(lfnlist.begin(),lfnlist.end(),DeleteLatticeFunction());
    lfnlist.clear();
    AddFunction(0,0,0);
    orbitonly = true;
};

struct ClearLatticeFunction
{
public:
    ClearLatticeFunction() {};
    void operator()(LatticeFunction* lfn) { lfn->ClearValues(); };
};

struct CalculateLatticeFunction
{
private:
    RealMatrix* N;
    const Particle* p;
    double s;

public:
    CalculateLatticeFunction(double _s, const Particle& _p, RealMatrix& _N)
            : s(_s), p(&_p), N(&_N) {};

    void operator()(LatticeFunction* lfn) {
        int i, j, k;
        double v = 0;

        lfn->GetIndices(i, j, k);

		if(i==0 && j==0 && k>0) {
			v = atan2( (*N)(2*k-2,2*k-1) , (*N)(2*k-2,2*k-2) )/twoPi;
			if(k!=3 && v<-1.0e-9)
				v += 1.0;
		}

        if(i==0 && j==0 && k==0)
            v = s;

        if(i>0 && j==0 && k==0)
            v = (*p)[i-1];

        if(i>0 && j>0 && k>0)
            v = (*N)(i-1,2*k-2)*(*N)(j-1,2*k-2) + (*N)(i-1,2*k-1)*(*N)(j-1,2*k-1);

        lfn->AppendValue(v);
    };
};

void LatticeFunctionTable::Calculate(PSvector* p, RealMatrix* M)
{
    if(orbitonly)
        DoCalculateOrbitOnly(bendscale, p);
    else
        DoCalculate(bendscale, p, M);
}

struct CopyLatticeFunction
{
private:
    vectorlfn* lfncopy;
public:
    CopyLatticeFunction(vectorlfn& copy) : lfncopy(&copy) {};

    void operator()(LatticeFunction* lfn){
        LatticeFunction* newlfn = new LatticeFunction(*lfn);
        lfncopy->push_back(newlfn);
    };
};

void LatticeFunctionTable::CalculateEnergyDerivative()
{
    double dpP = orbitonly ? DoCalculateOrbitOnly( bendscale) : DoCalculate( bendscale);
    vectorlfn lfnP;
    for_each(lfnlist.begin(),lfnlist.end(),CopyLatticeFunction(lfnP));

    double dpM = orbitonly ? DoCalculateOrbitOnly(-bendscale) : DoCalculate(-bendscale);
    vectorlfn lfnM;
    for_each(lfnlist.begin(),lfnlist.end(),CopyLatticeFunction(lfnM));

    for_each(lfnlist.begin(),lfnlist.end(),ClearLatticeFunction());

    double dp = dpP - dpM;
    vectorlfn::iterator lfnitP = lfnP.begin();
    vectorlfn::iterator lfnitM = lfnM.begin();

    for(vectorlfn::iterator lfnit=lfnlist.begin(); lfnit!=lfnlist.end(); lfnit++,lfnitP++,lfnitM++) {
        (*lfnit)->Derivative(*lfnitM, *lfnitP, dp);
    };

}

double LatticeFunctionTable::DoCalculate(double cscale, PSvector* pInit, RealMatrix* MInit)
{
    for_each(lfnlist.begin(),lfnlist.end(),ClearLatticeFunction());

    PSvector p(0);
	if(pInit) {
		p = *pInit;
	} else {
		ClosedOrbit co(theModel, p0);
		co.SetDelta(delta);
		co.ScaleBendPathLength(cscale);
		co.FindClosedOrbit(p);
	}

    RealMatrix M(6);
	if(MInit) {
		M = *MInit;
	} else {
		TransferMatrix tm(theModel, p0);
		tm.SetDelta(delta);
		tm.ScaleBendPathLength(cscale);
		tm.FindTM(M,p);
	}
	
    ComplexVector eigenvalues(3);
    ComplexMatrix eigenvectors(3,6);
	if(symplectify)
		Symplectify(M);
    EigenSystem(M, eigenvalues, eigenvectors);

	int row, col;

    RealMatrix N(6);
	RealMatrix R(6);
    for(row=0; row<6; row++) {
        for(col=0; col<3; col++) {
            N(row,2*col)   = sqrt(2.)*eigenvectors(col,row).real();
            N(row,2*col+1) = sqrt(2.)*eigenvectors(col,row).imag();
			R(row,2*col)   = 0.0;
			R(row,2*col+1) = 0.0;
        }
    }

	for(row=0; row<3; row++) {
		int i = 2*row;
		int j = 2*row + 1;
		double theta =-atan2(N(i,j),N(i,i));

		R(i,i) = R(j,j) = cos(theta);
		R(i,j) = sin(theta);
		R(j,i) =-sin(theta);
	}

	N = N*R;

//	ofstream nfile("DataFiles/NormMatrix.dat");
//	MatrixForm(N,nfile,OPFormat().precision(6).fixed());

    ParticleBunch* particle = new ParticleBunch(p0, 1.0);
    particle->push_back(p);
    particle->push_back(p);

    for(row=0; row<6; row++) {
        Particle q = p;
        q[row] += delta;
        particle->push_back(q);
    };

    ParticleTracker tracker(theModel->GetBeamline(),particle);

    if(cscale) {
        RingDeltaTProcess* ringdt = new RingDeltaTProcess(2);
        ringdt->SetBendScale(cscale);
        tracker.AddProcess(ringdt);
    }

    bool loop = true;
	bool isMore = true;
    tracker.InitStepper();

    RealMatrix M1 = IdentityMatrix(6);
	RealMatrix M2(6);
	RealMatrix M21(6);

	double e0 = particle->GetReferenceMomentum();
	double e1 = e0;
    double s  = 0;

	do {

		ParticleBunch::iterator ip = tracker.GetTrackedBunch().begin();
		Particle pref1 = *ip++;
		Particle pref2 = *ip++;

		e1 = tracker.GetTrackedBunch().GetReferenceMomentum();

		if(e1!=e0)
//			for(int row=0; row<3; row++)
//				pref2[2*row+1] *= e1/e0;
			for(int row=0; row<6; row++)
				pref2[row] *= sqrt(e1/e0);
			
		for(int col=0; col<6; col++,ip++) {

			if(e1!=e0)
//				for(int row=0; row<3; row++)
//					(*ip)[2*row+1] *= e1/e0;
				for(int row=0; row<6; row++)
					(*ip)[row] *= sqrt(e1/e0);

			for(int row=0; row<6; row++) {
				M2(row,col) = ((*ip)[row] - pref2[row]) / delta;
			}
		}

		e0 = e1;

		M21 = M2*M1;
		M1  = M2;
		Invert(M1);

		if(symplectify)
			Symplectify(M21);

		N  = M21*N;

		for_each(lfnlist.begin(),lfnlist.end(),CalculateLatticeFunction(s,pref1,N));
		if(isMore){
			s += tracker.GetCurrentComponent().GetLength();
			isMore = tracker.StepComponent();
		}else
			loop = false;

    } while(loop);

//	ofstream mfile("TransferMatrix.dat");
//	MatrixForm(M2,mfile,OPFormat().precision(6).fixed());

    return p.dp();
};

double LatticeFunctionTable::DoCalculateOrbitOnly(double cscale, PSvector* pInit)
{
    for_each(lfnlist.begin(),lfnlist.end(),ClearLatticeFunction());

    PSvector p(0);
	if(pInit) {
		p = *pInit;
	} else {
		ClosedOrbit co(theModel, p0);
		co.SetDelta(delta);
		co.ScaleBendPathLength(cscale);
		co.FindClosedOrbit(p);
	}

    ParticleBunch* particle = new ParticleBunch(p0, 1.0);
    particle->push_back(p);

    ParticleTracker tracker(theModel->GetBeamline(),particle);

    if(cscale) {
        RingDeltaTProcess* ringdt = new RingDeltaTProcess(2);
        ringdt->SetBendScale(cscale);
        tracker.AddProcess(ringdt);
    }

    bool loop = true;
    tracker.InitStepper();

    double s = 0;
    RealMatrix N1(6);

    do {

        ParticleBunch::const_iterator ip = tracker.GetTrackedBunch().begin();
        const Particle& pref = *ip++;

        for_each(lfnlist.begin(),lfnlist.end(),CalculateLatticeFunction(s,pref,N1));
        s += tracker.GetCurrentComponent().GetLength();
        loop = tracker.StepComponent();

    } while(loop);

    return p.dp();
};

struct PrintLatticeFunction
{
private:
    ostream* os;
    int row;

public:
    PrintLatticeFunction(int r, ostream& _os)
            : row(r), os(&_os) {};

    void operator()(LatticeFunction* lfn) {
        (*os)<<std::setw(14)<<lfn->GetValue(row);
    };

};

inline int LatticeFunctionTable::NumberOfRows()
{
    return (*(lfnlist.begin()))->size();
};

void LatticeFunctionTable::PrintTable(ostream& os, int n1, int n2)
{
    if(n1<0)
        n1 = 0;

    int nrows = NumberOfRows();
    if(n2>=nrows || n2<n1)
        n2 = nrows-1;

    for(int row=n1; row<=n2; row++) {
        for_each(lfnlist.begin(), lfnlist.end(), PrintLatticeFunction(row,os));
        os<<endl;
    };
};

int LatticeFunctionTable::GetSPosIndex(double s)
{
    int n;

    for(n=0; n<NumberOfRows(); n++)
        if(lfnlist[0]->GetValue(n)>=s)
            return n;

    return -1;
};

double LatticeFunctionTable::Value(int i, int j, int k, int ncpt)
{
    vectorlfn::iterator lfnit = GetColumn(i, j, k);
    return (*lfnit)->GetValue(ncpt);
};

vectorlfn::iterator LatticeFunctionTable::GetColumn(int i, int j, int k)
{
    vectorlfn::iterator lfnit = lfnlist.begin();
    int _i, _j, _k;

    for(;lfnit!=lfnlist.end();lfnit++) {
        (*lfnit)->GetIndices(_i,_j,_k);
        if(i==_i && j==_j && k==_k)
            return lfnit;
    };

    return lfnit;
};

double LatticeFunctionTable::Mean(int i, int j, int k, int n1, int n2)
{
    vectorlfn::iterator lfnit = GetColumn(i,j,k);

    if(n1<0)
        n1 = 0;

    int nrows = NumberOfRows();
    if(n2>=nrows || n2<n1)
        n2 = nrows-1;

    double sums = 0;
    double sumv = 0;
    for(int row=n1; row<n2; row++) {
        double ds = lfnlist[0]->GetValue(row+1) - lfnlist[0]->GetValue(row);
        double mv = (*lfnit)->GetValue(row+1) + (*lfnit)->GetValue(row);
        sums += ds;
        sumv += 0.5*ds*mv;
    };

    return sumv/sums;
};

double LatticeFunctionTable::RMS(int i, int j, int k, int n1, int n2)
{
    vectorlfn::iterator lfnit = GetColumn(i,j,k);

    if(n1<0)
        n1 = 0;

    int nrows = NumberOfRows();
    if(n2>=nrows || n2<n1)
        n2 = nrows-1;

    double sums = 0;
    double sumv = 0;
    for(int row=n1; row<n2; row++) {
        double ds = lfnlist[0]->GetValue(row+1) - lfnlist[0]->GetValue(row);
        double mv1 = (*lfnit)->GetValue(row+1);
        double mv0 = (*lfnit)->GetValue(row);
        sums += ds;
        sumv += 0.5*ds*(mv0*mv0 + mv1*mv1);
    };

    return sqrt(sumv/sums);
};
