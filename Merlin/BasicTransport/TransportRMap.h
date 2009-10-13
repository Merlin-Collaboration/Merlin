/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef TransportRMap_h
#define TransportRMap_h 1

#include "merlin_config.h"
#include "BasicTransport/RMap.h"

namespace TransportRMap {

void Drift (double length, RMap& R);
void SectorBend (double length, double h, double Kx, RMap& R);
void SectorBendT (double length, double h, double Kx, RMap& T);
void Quadrupole (double length, double Kx, RMap& R);
void Srot (double phi, RMap& R);
void Srot (double cosphi, double sinphi, RMap& R);
void PoleFaceRot (double h, double theta, double fint, double hgap, RMap& R);
void Solenoid (double length, double K0, double K1, bool entrFringe, bool exitFringe, RMap& R);
void TWRFCavity (double length, double g, double f, double phi,
                 double E0, bool entr_field, bool exit_field, RMap& R);
void SWRFCavity (int ncells, double g, double f, double phi, double E0, RMap& R);

inline void Srot (double phi, RMap& R) {
    Srot(cos(phi),sin(phi),R);
}

// functions returning R2Map objects
void Drift (double length, R2Map& R);
void Quadrupole (double length, double Kx, R2Map& R);
void TWRFCavity (double length, double g, double f, double phi,
                 double E0, bool entr_field, bool exit_field, R2Map& R);
void SWRFCavity (int ncells, double g, double f, double phi, double E0, R2Map& R);

}; // namespace TransportRMap

#endif
