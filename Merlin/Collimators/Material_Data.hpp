#ifndef _Material_Data_hpp_
#define _Material_Data_hpp_

// proton nucleus total cross section =sigma_el_pn+sigma_el_pN+sigma_SD_pn +sigma_inel_pN 
//                                           Be     C      Al     Cu     W      Pb
const double coll_mat_data_sig_pN_tot[] = {0.268, 0.331, 0.634, 1.232, 2.767, 2.960};

// proton nucleus inelastic cross section (sigma_inel_pN) 
//                                           Be     C      Al     Cu     W      Pb
const double coll_mat_data_sig_pN_in[] = {0.199, 0.231, 0.421, 0.782, 1.65, 1.77};

// Rutherford scattering cross section (sigma_R) 
//					Be		C		Al		Cu		W		Pb
const double coll_mat_data_sig_R[] = 	{0.000035,	0.000076,	0.00034,	0.00153,	0.00768,	0.00907};
const double coll_mat_data_dEdx[] = 	{1.594,		1.745,		1.724,		1.403,		1.145,		1.123};
const double coll_mat_data_rho[] = 	{1.848,		2.265,		1.204,		8.96,		19.3,		11.35}; 
const double coll_mat_data_A[] = 	{9.012182,	12.011,		1.204,		26.981539,	183.84,		207.2}; 

// electrical conductivity (sigma) =1/electrical resisitivity (Ohm*m)e-1
//                                     Be        C     Al       Cu       W        Pb       
const double coll_mat_data_sigma[] = {3.08e7, 7.14e4, 35.64e6, 5.98e7, 0.177e4, 4.8077e6 };
// real-space radiation length (m)
//					Be        C     Al       Cu       W        Pb 
const double coll_mat_data_X0[] = 	{35.28e-2, 18.8e-2, 8.90e-2, 1.44e-2, 0.351e-2, 0.562e-2 };

extern const char* material_names[];

#endif
