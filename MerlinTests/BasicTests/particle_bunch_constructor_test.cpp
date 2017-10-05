#include <iostream>
#include <tuple>

#include "../tests.h"
#include "ParticleBunchTypes.h"
#include "ParticleBunchConstructor.h"
#include "BeamData.h"
#include "ParticleBunch.h"
#include "RandomNG.h"

using namespace std;

/*
 * Construct bunches with a range of distributions and beam parameters
 *
 * Check that the bunch is the correct size
 * Check the initial particle is on the centroid
 * Check the mean is close to the centroid
 * Check the standard deviations for flat and normal
 * Check the extents for flat and ring
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

	for (auto & p: *pb)
	{
		for(int i=0; i<PS_LENGTH ; i++)
		{
			bunch_stats.min[i] = min(p[i], bunch_stats.min[i]);
			bunch_stats.max[i] = max(p[i], bunch_stats.max[i]);
			sum[i] += p[i];
			sum_sqs[i] += p[i]*p[i];
		}
	}

	int npart = pb->GetParticles().size();
	for(int i=0; i<PS_LENGTH ; i++)
	{
		bunch_stats.mean[i] = sum[i]/ npart;
		bunch_stats.std[i] = sqrt((sum_sqs[i]  - sum[i]*sum[i]/npart)/ npart);
	}

	return bunch_stats;
}

// compare PSvector to values, with tolerance.
bool are_close(PSvector a, double x, double xp, double y, double yp, double ct, double dp, double tol)
{
	bool good = true;
	if (fabs(a.x() - x) > tol)
	{
		cout << "are_close() x: " << a.x() << "!=" << x << endl;
		good = false;
	}
	if (fabs(a.xp() - xp) > tol)
	{
		cout << "are_close() xp: " << a.xp() << "!=" << xp << endl;
		good = false;
	}
	if (fabs(a.y() - y) > tol)
	{
		cout << "are_close() y: " << a.y() << "!=" << y << endl;
		good = false;
	}
	if (fabs(a.yp() - yp) > tol)
	{
		cout << "are_close() yp: " << a.yp() << "!=" << yp << endl;
		good = false;
	}
	if (fabs(a.ct() - ct) > tol)
	{
		cout << "are_close() ct: " << a.ct() << "!=" << ct << endl;
		good = false;
	}
	if (fabs(a.dp() - dp) > tol)
	{
		cout << "are_close() dp: " << a.dp() << "!=" << dp << endl;
		good = false;
	}
	return good;
}

int main()
{
	size_t npart = (size_t) 1e6;

	int seed = 1;
	RandomNG::init(seed);

	auto dists = {make_tuple("normal", normalDistribution),
	              make_tuple("flat", flatDistribution),
	              make_tuple("ring", ringDistribution),
	              make_tuple("horizontalHalo", horizontalHaloDistribution1),
	              make_tuple("verticalHalo", verticalHaloDistribution1),
	              make_tuple("horizontalHalo2", horizontalHaloDistribution2),
	              make_tuple("verticalHalo2", verticalHaloDistribution2)
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

	for(auto & beam: beams)
	{
		cout << endl << "x0: " << beam.x0
		     << " xp0: " << beam.xp0
		     << " y0: " << beam.y0
		     << " yp0: " << beam.yp0
		     << " emit_x: " << beam.emit_x
		     << " emit_y: " << beam.emit_y
		     << " beta_x:" << beam.beta_x
		     << " beta_y: " << beam.beta_y
		     << endl;

		for (auto& dist :dists)
		{
			auto name = get<0>(dist);
			auto dist_type = get<1>(dist);

			cout << "Dist: " << name << endl;
			ParticleBunchConstructor pbc(beam, npart, dist_type);
			ParticleBunch* myBunch = pbc.ConstructParticleBunch<ParticleBunch>();

			assert(myBunch->GetParticles().size() == npart);

			// Check first particle is centroid
			auto &p0 = myBunch->GetParticles()[0];
			assert(are_close(p0, beam.x0, beam.xp0, beam.y0, beam.yp0, 0, 0, 0));

			stats bunch_stats = get_stats(myBunch);

			// check mean
			assert(are_close(bunch_stats.mean, beam.x0, beam.xp0, beam.y0, beam.yp0, 0, 0, 1e-2));

			// check standard deviations
			double gamma_x = (1 + beam.alpha_x*beam.alpha_x) / beam.beta_x;
			double gamma_y = (1 + beam.alpha_y*beam.alpha_y) / beam.beta_y;
			double sig_x = sqrt(beam.emit_x * beam.beta_x);
			double sig_xp = sqrt(beam.emit_x * gamma_x);
			double sig_y = sqrt(beam.emit_y * beam.beta_y);
			double sig_yp = sqrt(beam.emit_y * gamma_y);

			if (name == string("flat"))
			{
				double c = sqrt(1.0/12.0)*2; // std of a uniform dist
				assert(are_close(bunch_stats.std, c*sig_x, c*sig_xp, c*sig_y, c*sig_yp, 0, 0, 1e-2));
			}
			if (name == string("normal"))
			{
				assert(are_close(bunch_stats.std, sig_x, sig_xp, sig_y, sig_yp, 0, 0, 1e-2));
			}

			// check extents
			if ((name == string("flat")) || (name == string("ring")))
			{
				assert(are_close(bunch_stats.max, beam.x0+sig_x, beam.xp0+sig_xp, beam.y0+sig_y, beam.yp0+sig_yp, 0, 0, 1e-3));
				assert(are_close(bunch_stats.min, beam.x0-sig_x, beam.xp0-sig_xp, beam.y0-sig_y, beam.yp0-sig_yp, 0, 0, 1e-3));
			}

			delete myBunch;
		}
	}
	cout << "Done" << endl;
	return 0;
}