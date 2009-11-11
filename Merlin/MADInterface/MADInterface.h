/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.7 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef MADInterface_h
#define MADInterface_h 1

#include "merlin_config.h"
#include <fstream>
#include <string>
#include <set>
// AcceleratorModel
#include "AcceleratorModel/AcceleratorModel.h"
//#include "AcceleratorModel/Components.h"
class AcceleratorModelConstructor;
class MADKeyMap;

using std::ifstream;
using std::ostream;

//      Class used to construct a MERLIN model from a MAD optics
//      output listing. The class now automatically  identifies
//      the column parameters, and associates them with the
//      constructed element types. If an element type is defined
//      for which a required parameter is not present in the
//      column headings, the parameter is set to zero and a
//      warning is issued.

class MADInterface
{
public:

    //  Constructor taking the name of the MAD optics file, and
    //  the momentum in GeV/c.
    MADInterface (const std::string& madFileName="", double P0=0);

    //   Causes the construction of an AcceleratorModel object
    //   based on the MAD optics file.
    AcceleratorModel* ConstructModel ();

    //   Sets the log file stream to os.
    void SetLogFile (ostream& os);

    //   Turns logging on.
    void SetLoggingOn ();

    //   Turns logging off.
    void SetLoggingOff ();

    //   If true, all LINE constructs in the MAD optics output
    //   are constructed in the model. If false, only those
    //   prefixed X_, where X is M, S, or G are constructed.
    void HonourMadStructure (bool flg);

    //   If true, a flat lattice model in constructed, with no
    //   nested frames.
    void ConstructFlatLattice (bool flg);

    void ConstructApertures (bool flg);

    //   Components of type madType are ignored during
    //   construction if their length is zero.
    void IgnoreZeroLengthType (const string& madType);

    //   If scaleSR == true, then the magnetic fields of the
    //   magnets are scaled to compensate beam energy losses due
    //   to synchrotron radiation (default = false.) Note that in
    //   this case, the beam energy is the initial energy.
    void ScaleForSynchRad (bool scaleSR);

    //   Treats the mad type typestr as a drift.
    void TreatTypeAsDrift (const std::string& typestr);

    // Functions for constructing a model from several files.
    // Repeated calls to AppendModel(fname,p) constructs a single
    // model (beamline) from the respective files. The final
    // model is returned using GetModel().
    void AppendModel(const std::string& fname, double pref);
    AcceleratorModel* GetModel();

    void ConstructNewFrame (const string& name);
    void EndFrame (const string& name);

protected:

    double energy;
    ifstream *ifs;
    ostream* log;

    bool logFlag;
    bool flatLattice;
    bool honMadStructs;
    bool incApertures;
    bool inc_sr;

    std::set<std::string> zeroLengths;
    std::set<std::string> driftTypes;

    AcceleratorModelConstructor* ctor;
    MADKeyMap* prmMap;

    double ReadComponent ();
    void Initialise();
};

inline void MADInterface::SetLogFile (ostream& os)
{
    log=&os;
}

inline void MADInterface::SetLoggingOn ()
{
    logFlag=true;
}

inline void MADInterface::SetLoggingOff ()
{
    logFlag=false;
}

inline void MADInterface::HonourMadStructure (bool flg)
{
    honMadStructs=flg;
}

inline void MADInterface::ConstructFlatLattice (bool flg)
{
    flatLattice=flg;
}

inline void MADInterface::ConstructApertures (bool flg)
{
    incApertures = flg;
}

inline void MADInterface::ScaleForSynchRad (bool scaleSR)
{
    inc_sr = scaleSR;
}



#include "AcceleratorModel/Apertures/SimpleApertures.h"

// proton nucleus total cross section =sigma_el_pn+sigma_el_pN+sigma_SD_pn +sigma_inel_pN 
//                                           Be     C      Al     Cu     W      Pb
const double coll_mat_data_sig_pN_tot[] = {0.268, 0.331, 0.634, 1.232, 2.767, 2.960};

// proton nucleus inelastic cross section (sigma_inel_pN) 
//                                           Be     C      Al     Cu     W      Pb
const double coll_mat_data_sig_pN_in[] = {0.199, 0.231, 0.421, 0.782, 1.65, 1.77};

// Rutherford scattering cross section (sigma_R) 
//                                        Be        C        Al       Cu     W        Pb
const double coll_mat_data_sig_R[] = {0.000035, 0.000076, 0.00034, 0.00153, 0.00768, 0.00907};
const double coll_mat_data_dEdx[] = {1.594, 1.745, 1.724, 1.403, 1.145, 1.123};
const double coll_mat_data_rho[] = {1.848, 2.265, 1.204, 8.96, 19.3, 11.35}; 
const double coll_mat_data_A[] = {9.012182, 12.011, 1.204, 26.981539, 183.84, 207.2}; 

// electrical conductivity (sigma) =1/electrical resisitivity (Ohm*m)e-1
//                                     Be        C     Al       Cu       W        Pb       
const double coll_mat_data_sigma[] = {3.08e7, 7.14e4, 35.64e6, 5.98e7, 0.177e4, 4.8077e6 };
// real-space radiation length (m)
//                                     Be        C     Al       Cu       W        Pb 
const double coll_mat_data_X0[] = {35.28e-2, 18.8e-2, 8.90e-2, 1.44e-2, 0.351e-2, 0.562e-2 };

 
extern const char* material_names[]; 


class TiltedAperture: public RectangularAperture
{ double alpha;
int material;
//xoff,yoff;

public:
double sig_pN_tot;
double sig_pN_in;
double sig_R;
double dEdx;
double rho;
double sigma;
double A;
double X0;
  //double tot_mean_free_path;
         TiltedAperture(double w,double h, double t, int m=0):RectangularAperture(w,h), alpha(t), material(m){

	   //  cout << "Tilted aperture made of " << material_names[material]<<endl;
		sig_pN_tot = coll_mat_data_sig_pN_tot[material];
		sig_pN_in = coll_mat_data_sig_pN_in[material];
		sig_R = coll_mat_data_sig_R[material];
                dEdx = coll_mat_data_dEdx[material];
                rho = coll_mat_data_rho[material];
                A = coll_mat_data_A[material];
	        X0=coll_mat_data_X0[material];

                sigma = coll_mat_data_sigma[material];
                
		//tot_mean_free_path = 1 / sig_pN_tot;
	};
         virtual bool PointInside(double x,double y,double z) const
         {
		    double x1 = x * cos(alpha);
		    double y1 = y * sin(alpha);
		 //   cout << "pointinside (" << x1 <<"," << y1 << ")"<<endl;
		    return fabs(x1)<GetFullWidth()/2 && fabs(y1)<GetFullHeight()/2;

	    };
       /*  virtual void Setoff(const double x,const double y)
         {
              xoff=x;
              yoff=y;
		}*/

		virtual bool PointInside_offset(double x,double y,double z,const double xoff,const double yoff) const
		         {
				    double x1 = (x-xoff) * cos(alpha);
				    double y1 = (y-yoff) * sin(alpha);
				    // cout << "pointinside (" << x1 <<"," << y1 << ")"<<endl;
	//			     cout << "offset (" << xoff <<"," << yoff << ")"<<endl;
				    return fabs(x1)<GetFullWidth()/2 && fabs(y1)<GetFullHeight()/2;

	    };
	};


#endif


