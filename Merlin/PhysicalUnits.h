/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_PhysicalUnits
#define _h_PhysicalUnits

// Merlin common units

namespace PhysicalUnits
{

//Area
extern const double barn;

// length
extern const double meter;
extern const double centimeter;
extern const double millimeter;
extern const double micrometer;
extern const double nanometer;

// time
extern const double second;
extern const double millisecond;
extern const double microsecond;
extern const double nanosecond;
extern const double picosecond;
extern const double minute;
extern const double hour;
extern const double day;
extern const double year;

// energy
extern const double eV;
extern const double keV;
extern const double MeV;
extern const double GeV;
extern const double TeV;

// voltage
extern const double Volt;
extern const double kV;
extern const double MV;
extern const double GV;
extern const double TV;

// Frequency
extern const double Hz;
extern const double kHz;
extern const double MHz;
extern const double GHz;
extern const double THz;

// Angle
extern const double radian;
extern const double milliradian;
extern const double microradian;
extern const double degree;

// magnetic field
extern const double Tesla;
extern const double kGauss;
extern const double Gauss;

// mass
extern const double amu;

}

#endif
