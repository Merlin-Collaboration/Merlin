/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "PSvectorTransform3D.h"

PSvectorTransform3D::PSvectorTransform3D(const Transform3D& tfrm) :
	T(tfrm), bNoRot(tfrm.R().isIdentity())
{
}

PSvector& PSvectorTransform3D::Apply(PSvector& p) const
{
	if(bNoRot)
	{
		// drift space with transverse displacement
		p.x() += (T.X().z * p.xp() - T.X().x);
		p.y() += (T.X().z * p.yp() - T.X().y);
		//	p.ct()+=T.X().z; //dubious
	}
	else
	{
		Point3D X(p.x(), p.y(), 0);
		Vector3D V(p.xp(), p.yp(), 1.0);

		X = T(X);
		V = T(V);

		// re-scale the "momentum" vector to (x',y',1)
		V.x /= V.z;
		V.y /= V.z;

		// copy back to p, applying any necessary drift
		// to bring z to zero
		p.x() = X.x - X.z * V.x;
		p.y() = X.y - X.z * V.y;
		//	p.ct()-= X.z; //dubious
		p.xp() = V.x;
		p.yp() = V.y;
	}
	return p;
}

PSvectorArray& PSvectorTransform3D::Apply(PSvectorArray& pv) const
{
	std::for_each(pv.begin(), pv.end(), *this);
	return pv;
}
