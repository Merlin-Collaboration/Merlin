/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <utility>
#include <fstream>

#include "../tests.h"
#include "ParticleBunchTypes.h"
#include "BeamData.h"
#include "ParticleBunch.h"
#include "RandomNG.h"
#include "ParticleDistributionGenerator.h"
#include "HaloParticleDistributionGenerator.h"

using namespace std;

/*
 * Construct bunches with a range of distributions and beam parameters
 *
 * Check that the bunch is the correct size
 * Check the initial particle is on the centroid
 * Check the mean is close to the centroid
 * Check the standard deviations for flat and normal
 * Check the extents for flat and ring
 * Check the ellipse for halo distributions
 *
 */

struct stats
{
	PSvector mean;
	PSvector std;
	PSvector min;
	PSvector max;

};

// Calculate min, max, mean and std deviation
stats get_stats(ParticleBunch * pb)
{
	stats bunch_stats;
	bunch_stats.min = pb->GetParticles()[0];
	bunch_stats.max = pb->GetParticles()[0];
	bunch_stats.mean = PSvector(0);
	PSvector sum(0);
	PSvector sum_sqs(0);

	for(auto & p : *pb)
	{
		for(int i = 0; i < PS_LENGTH; i++)
		{
			bunch_stats.min[i] = min(p[i], bunch_stats.min[i]);
			bunch_stats.max[i] = max(p[i], bunch_stats.max[i]);
			sum[i] += p[i];
			sum_sqs[i] += p[i] * p[i];
		}
	}

	int npart = pb->GetParticles().size();
	for(int i = 0; i < PS_LENGTH; i++)
	{
		bunch_stats.mean[i] = sum[i] / npart;
		bunch_stats.std[i] = sqrt((sum_sqs[i]  - sum[i] * sum[i] / npart) / npart);
	}

	return bunch_stats;
}

bool are_close(double a, double b, double tol)
{
	return fabs(a - b) < tol;
}

// compare PSvector to values, with tolerance.
bool are_close(PSvector a, double x, double xp, double y, double yp, double ct, double dp, double tol)
{
	bool good = true;
	if(fabs(a.x() - x) > tol)
	{
		cout << "are_close() x: " << a.x() << "!=" << x << endl;
		good = false;
	}
	if(fabs(a.xp() - xp) > tol)
	{
		cout << "are_close() xp: " << a.xp() << "!=" << xp << endl;
		good = false;
	}
	if(fabs(a.y() - y) > tol)
	{
		cout << "are_close() y: " << a.y() << "!=" << y << endl;
		good = false;
	}
	if(fabs(a.yp() - yp) > tol)
	{
		cout << "are_close() yp: " << a.yp() << "!=" << yp << endl;
		good = false;
	}
	if(fabs(a.ct() - ct) > tol)
	{
		cout << "are_close() ct: " << a.ct() << "!=" << ct << endl;
		good = false;
	}
	if(fabs(a.dp() - dp) > tol)
	{
		cout << "are_close() dp: " << a.dp() << "!=" << dp << endl;
		good = false;
	}
	return good;
}

int main(int argc, char* argv[])
{
	size_t npart = (size_t) 1e6;

	int seed = 1;
	RandomNG::init(seed);

	string ref_file_name = "";
	ofstream *ref_file;
	for(int i = 1; i < argc - 1; i++)
	{
		if(strcmp(argv[i], "--reffile") == 0)
		{
			ref_file_name = argv[i + 1];
		}
	}

	if(ref_file_name != "")
	{
		ref_file = new ofstream(ref_file_name);
		if(!ref_file->good())
		{
			cerr << "Failed to open reference file:" << ref_file_name << endl;
		}
		cout << "Writing reference file: " << ref_file_name << endl;
		ref_file->setf(ios::scientific, ios::floatfield);
		ref_file->precision(16);
		npart = 100;
	}

	double cut = 1.3;

	vector<pair<string, ParticleDistributionGenerator*> > dists =
	{
		{"normal", new NormalParticleDistributionGenerator()},
		{"normal_cut", new NormalParticleDistributionGenerator(cut)},
		{"normal_cent", new NormalParticleDistributionGenerator()},
		{"flat", new UniformParticleDistributionGenerator()},
		{"ring", new RingParticleDistributionGenerator()},
		{"horizontalHalo", new HorizonalHalo1ParticleDistributionGenerator()},
		{"verticalHalo", new VerticalHalo1ParticleDistributionGenerator()},
		{"horizontalHalo2", new HorizonalHalo2ParticleDistributionGenerator()},
		{"verticalHalo2", new VerticalHalo2ParticleDistributionGenerator()},
		{"skewHalo", new RingParticleDistributionGenerator()}, // note ring and skewHalo are identical in old code
	};

	vector<BeamData> beams;
	BeamData basebeam;
	basebeam.emit_x = 1;
	basebeam.emit_y = 1;
	basebeam.beta_x = 1;
	basebeam.beta_y = 1;
	basebeam.p0 = 1;
	beams.push_back(basebeam);

	BeamData abeam(basebeam);
	abeam.emit_x = 2;
	abeam.emit_y = 3;
	beams.push_back(abeam);

	abeam = basebeam;
	abeam.beta_x = 4;
	abeam.beta_y = 3;
	beams.push_back(abeam);

	abeam = basebeam;
	abeam.x0 = -2;
	abeam.y0 = 3;
	beams.push_back(abeam);

	abeam = basebeam;
	abeam.xp0 = -2;
	abeam.yp0 = 3;
	beams.push_back(abeam);

	abeam = basebeam;
	abeam.emit_x = 2;
	abeam.emit_y = 3;
	abeam.xp0 = -2;
	abeam.yp0 = 3;
	beams.push_back(abeam);

	for(auto & beam : beams)
	{
		cout << endl << "x0: " << beam.x0
			 << " xp0: " << beam.xp0
			 << " y0: " << beam.y0
			 << " yp0: " << beam.yp0
			 << " emit_x: " << beam.emit_x
			 << " emit_y: " << beam.emit_y
			 << " beta_x:" << beam.beta_x
			 << " beta_y: " << beam.beta_y
			 << " alpha_x:" << beam.alpha_x
			 << " alpha_y: " << beam.alpha_y
			 << " gamma_x:" << beam.gamma_x()
			 << " gamma_y: " << beam.gamma_y()
			 << endl;

		if(ref_file_name != "")
		{
			*ref_file << endl << "#x0: " << beam.x0
					  << " xp0: " << beam.xp0
					  << " y0: " << beam.y0
					  << " yp0: " << beam.yp0
					  << " emit_x: " << beam.emit_x
					  << " emit_y: " << beam.emit_y
					  << " beta_x:" << beam.beta_x
					  << " beta_y: " << beam.beta_y
					  << endl;
		}

		for(auto& dist :dists)
		{
			auto name = get<0>(dist);
			auto dist_type = get<1>(dist);

			cout << "Dist: " << name << endl;

			ParticleBunch* myBunch = new ParticleBunch(npart, *dist_type, beam);

			if(name == string("normal_cent"))
			{
				myBunch->SetCentroid();
			}

			if(ref_file_name != "")
			{
				*ref_file << "# Dist: " << name << endl;
				*ref_file << "#X XP Y YP CT DP TYPE LOCATION ID SD" << endl;
				for(auto &p : *myBunch)
				{
					for(int i = 0; i < PS_LENGTH; i++)
					{
						*ref_file << std::setw(35) << p[i];
					}
					*ref_file << endl;
				}
				continue;
			}

			assert(myBunch->GetParticles().size() == npart);

			// Check first particle is centroid
			auto &p0 = myBunch->GetParticles()[0];
			assert(are_close(p0, beam.x0, beam.xp0, beam.y0, beam.yp0, 0, 0, 0));

			stats bunch_stats = get_stats(myBunch);

			// check mean
			assert(are_close(bunch_stats.mean, beam.x0, beam.xp0, beam.y0, beam.yp0, 0, 0, 1e-2));

			if(name == string("normal_cent"))
			{
				assert(are_close(bunch_stats.mean, beam.x0, beam.xp0, beam.y0, beam.yp0, 0, 0, 1e-8));
			}

			// check standard deviations
			double gamma_x = (1 + beam.alpha_x * beam.alpha_x) / beam.beta_x;
			double gamma_y = (1 + beam.alpha_y * beam.alpha_y) / beam.beta_y;
			double sig_x = sqrt(beam.emit_x * beam.beta_x);
			double sig_xp = sqrt(beam.emit_x * gamma_x);
			double sig_y = sqrt(beam.emit_y * beam.beta_y);
			double sig_yp = sqrt(beam.emit_y * gamma_y);

			if(name == string("flat"))
			{
				double c = sqrt(1.0 / 12.0) * 2; // std of a uniform dist
				assert(are_close(bunch_stats.std, c * sig_x, c * sig_xp, c * sig_y, c * sig_yp, 0, 0, 1e-2));
			}
			if(name == string("normal"))
			{
				assert(are_close(bunch_stats.std, sig_x, sig_xp, sig_y, sig_yp, 0, 0, 1e-2));
			}

			// check extents
			if((name == string("flat")) || (name == string("ring")))
			{
				assert(are_close(bunch_stats.max, beam.x0 + sig_x, beam.xp0 + sig_xp, beam.y0 + sig_y, beam.yp0
					+ sig_yp, 0, 0, 1e-3));
				assert(are_close(bunch_stats.min, beam.x0 - sig_x, beam.xp0 - sig_xp, beam.y0 - sig_y, beam.yp0
					- sig_yp, 0, 0, 1e-3));
			}

			if(name == string("normal_cut"))
			{
				assert(are_close(bunch_stats.max, beam.x0 + sig_x * cut, beam.xp0 + sig_xp * cut, beam.y0 + sig_y * cut,
					beam.yp0 + sig_yp * cut, 0, 0, 1e-3));
				assert(are_close(bunch_stats.min, beam.x0 - sig_x * cut, beam.xp0 - sig_xp * cut, beam.y0 - sig_y * cut,
					beam.yp0 - sig_yp * cut, 0, 0, 1e-3));
			}

			if(name == string("horizontalHalo") || name == string("horizontalHalo2") || name == string("skewHalo"))
			{
				for(auto p = myBunch->begin() + 1; p != myBunch->end(); ++p)
				{
					double x = p->x() - beam.x0, xp = p->xp() - beam.xp0;
					double rx2 = x * x * gamma_x + xp * xp * beam.beta_x + 2 * x * xp * beam.alpha_x;
					assert(are_close(rx2, beam.emit_x, 1e-8));
				}
			}

			if(name == string("verticalHalo") || name == string("verticalHalo2") || name == string("skewHalo"))
			{
				for(auto p = myBunch->begin() + 1; p != myBunch->end(); ++p)
				{
					double y = p->y() - beam.y0, yp = p->yp() - beam.yp0;
					double ry2 = y * y * gamma_y + yp * yp * beam.beta_y + 2 * y * yp * beam.alpha_y;
					assert(are_close(ry2, beam.emit_y, 1e-8));
				}
			}

			delete myBunch;
		}
	}

	for(auto& dist :dists)
	{
		delete get<1>(dist);
	}
	cout << "Done" << endl;
	return 0;
}
