/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _collimatortable_h_
#define _collimatortable_h_

#include <iostream>
#include <fstream>
using namespace std;

class collimatortable
{
	double* coeff;
	double step; // step size in s
	int ncoeff; // number
	double lo, hi;
public:

	~collimatortable()
	{
		delete[] coeff;
	}

	bool inrange(double x)
	{
		return (x >= lo) && (x <= hi);
	}

	double parabolic(double v1, double v2, double v3, double d)
	{
		return v2 + d * (v3 - v1) / 2 + d * d * (v3 + v1 - 2 * v2) / 2;
	}

	double interpolate(double s)
	{
		int index = int(0.5 + s / step);

		if(index < 0)
		{
			index = 0;
		}
		int i1;

		if(index == 0)
		{
			i1 = 1;
		}
		else if(index == ncoeff - 1)
		{
			i1 = ncoeff - 2;
		}
		else
		{
			i1 = index;
		}

		return (index >= ncoeff) ? 0 : parabolic(coeff[i1 - 1], coeff[i1], coeff[i1 + 1], s / step - i1);
	}

	collimatortable(const char* file, double Gamma = 0, double xi = 0)
	{
		double lo1, hi1, lo3, hi3;
		ifstream f;
		f.open(file);
		if(!f)
		{
			cout << "Cannot open file " << file << endl;
			exit(1);
		}

		int n1, n2, n3;
		f >> n3 >> lo3 >> hi3;
		double step3 = (hi3 - lo3) / (n3 - 1);
		f >> n2 >> lo >> hi;
		step = (hi - lo) / (n2 - 1);
		f >> n1 >> lo1 >> hi1;
		double step1 = (hi1 - lo1) / (n1 - 1);
		ncoeff = n2;
		coeff = new double[n2];

		double ***array;
		array = new double**[n1];
		for(int ii = 0; ii < n1; ii++)
		{
			array[ii] = new double*[n2];
			for(int jj = 0; jj < n2; jj++)
			{
				array[ii][jj] = new double[n3];
			}
		}

		while(f.get() != '{')
			;
		for(int i1 = 0; i1 < n1; i1++)
		{
			while(f.get() != '{')
				;
			for(int i2 = 0; i2 < n2; i2++)
			{
				while(f.get() != '{')
					;
				for(int i3 = 0; i3 < n3; i3++)
				{
					f >> array[i1][i2][i3];
					if(i3 < n3 - 1)
						while(f.get() != ',')
							;
				}
				while(f.get() != '}')
					;
				if(i2 < n2 - 1)
					while(f.get() != ',')
						;
			}
			while(f.get() != '}')
				;
			if(i1 < n1 - 1)
				while(f.get() != ',')
					;
		}
		while(f.get() != '}')
			;

		if((Gamma == 0) && (xi == 0))   // In principle the interpolation should
		{
			// handle this but it seems cleaner this way
			for(int i = 0; i < ncoeff; i++)
			{
				coeff[i] = array[0][i][0];
			}
		}
		else if(Gamma == 0) // interpolate xi value
		{
			int index = int(0.5 + xi / step1);
			if(index < 0)
			{
				index = 0;
			}
			int i1;
			if(index == 0)
			{
				i1 = 1;
			}
			else if(index == n1 - 1)
			{
				i1 = n1 - 2;
			}
			else
			{
				i1 = index;
			}
			for(int i = 0; i < ncoeff; i++)
			{
				if(index >= n1)
				{
					coeff[i] = 0;
				}
				else
				{
					double d = (xi - i1 * step1) / step1;
					coeff[i] = parabolic(array[0][i][i1 - 1], array[0][i][i1], array[0][i][i1 + 1], d);
				}
			}
		}
		else if(xi == 0) // interpolate Gamma value
		{
			int index = int(0.5 + Gamma / step3);
			if(index < 0)
			{
				index = 0;
			}
			int i1;

			if(index == 0)
			{
				i1 = 1;
			}
			else if(index == n3 - 1)
			{
				i1 = n3 - 2;
			}
			else
			{
				i1 = index;
			}

			for(int i = 0; i < ncoeff; i++)
			{
				if(index >= n1)
				{
					coeff[i] = 0;
				}
				else
				{
					double d = (Gamma - i1 * step3) / step3;
					coeff[i] = parabolic(array[i1 - 1][i][0], array[i1][i][0], array[i1 + 1][i][0], d);
				}
			}
		}
		else // 2D interpolation
		{
			int indg = int(0.5 + Gamma / step3);
			if(indg < 0)
			{
				indg = 0;
			}
			int ig1;
			if(indg == 0)
			{
				ig1 = 1;
			}
			else if(indg == n1 - 1)
			{
				ig1 = n1 - 2;
			}
			else
			{
				ig1 = indg;
			}
			int indc = int(0.5 + xi / step1);
			if(indc < 0)
			{
				indc = 0;
			}
			int ic1;
			if(indc == 0)
			{
				ic1 = 1;
			}
			else if(indc == n3 - 1)
			{
				ic1 = n3 - 2;
			}
			else
			{
				ic1 = indc;
			}
			double dg = (Gamma - ig1 * step3) / step3;
			double dc = (xi - ic1 * step1) / step1;
			for(int i = 0; i < ncoeff; i++)
			{
				double f0 = array[ig1][i][ic1];
				double fc = (array[ig1][i][ic1 + 1] - array[ig1][i][ic1 - 1]) / 2;
				double fg = (array[ig1 + 1][i][ic1] - array[ig1 - 1][i][ic1]) / 2;
				double fcc = (array[ig1][i][ic1 + 1] + array[ig1][i][ic1 - 1]) / 2 - array[ig1][i][ic1];
				double fgg = (array[ig1 + 1][i][ic1] + array[ig1 - 1][i][ic1]) / 2 - array[ig1][i][ic1];
				double apm = array[ig1 + 1][i][ic1 - 1] - (f0 + fg + fgg - fc + fcc);
				double app = array[ig1 + 1][i][ic1 + 1] - (f0 + fg + fgg + fc + fcc);
				double amm = array[ig1 - 1][i][ic1 - 1] - (f0 - fg + fgg - fc + fcc);
				double amp = array[ig1 - 1][i][ic1 + 1] - (f0 - fgg + fc + fcc);
				double fgc = (app + amm - apm - amp) / 4.;
				double fggc = (app - amp + apm - amm) / 4.;
				double fgcc = (app + amp - apm - amm) / 4.;
				double fggcc = (apm + app + amp + amm) / 4.;

				coeff[i] = f0 + dg * fg + dc * fc + dg * dg * fgg + dc * dc * fcc + dg * dc * fgc + dg * dg * dc * fggc
					+ dg * dc * dc * fgcc + dg * dg * dc * dc * fggcc;
			}
		}
	}
};

#endif
