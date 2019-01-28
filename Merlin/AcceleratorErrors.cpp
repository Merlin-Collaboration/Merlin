/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "AcceleratorErrors.h"
#include "utils.h"

using namespace std;
namespace
{

class Errors
{
public:

	Errors(double vvx, double vvy, double vvz, double meanx, double meany, double meanz, const string& p, bool clear,
		bool trans, ostream* l) :
		vx(vvx), vy(vvy), vz(vvz), mx(meanx), my(meany), mz(meanz), c(clear), t(trans), pat("*." + p), log(l)
	{
	}

	void operator()(LatticeFrame* frame) const
	{
		if(frame && pat((*frame).GetQualifiedName()))
		{

			if(c)
			{
				frame->ClearLocalFrameTransform();
			}

			double ex = fequal(vx, 0.0) ? mx : RandomNG::normal(mx, vx);
			double ey = fequal(vy, 0.0) ? my : RandomNG::normal(my, vy);
			double ez = fequal(vz, 0.0) ? mz : RandomNG::normal(mz, vz);
			if(t)
			{
				frame->Translate(ex, ey, ez);

				if(log)
				{
					(*log) << (*frame).GetQualifiedName() << " translate: " << ex << " " << ey << " " << ez << endl;
				}
			}
			else
			{
				if(!fequal(ex, 0.0))
				{
					frame->RotateX(ex);
				}
				if(!fequal(ey, 0.0))
				{
					frame->RotateY(ey);
				}
				if(!fequal(ez, 0.0))
				{
					frame->RotateZ(ez);
				}

				if(log)
				{
					(*log) << (*frame).GetQualifiedName() << " rotate: " << ex << " " << ey << " " << ez << endl;
				}
			}
		}
	}

private:

	double vx, vy, vz;
	double mx, my, mz;
	bool c, t;
	StringPattern pat;
	ostream* log;
};
} // End namespace

void AcceleratorErrors::ApplyShifts(AcceleratorModel::Beamline& b, const string& p)
{
	for_each(b.begin(), b.end(), Errors(vx, vy, vz, mx, my, mz, p, clear, true, log));
}

void AcceleratorErrors::ApplyRotations(AcceleratorModel::Beamline& b, const string& p)
{
	for_each(b.begin(), b.end(), Errors(vx, vy, vz, mx, my, mz, p, clear, false, log));
}
