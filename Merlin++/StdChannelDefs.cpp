/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "LatticeFrame.h"
#include "MagnetMover.h"
#include "BPMChannelCtor.h"
#include "CorrectorWinding.h"
#include "Components.h"
#include "IndirectChannels.h"

#define DEF_RW_CHANNEL(type, key, getf, setf) \
	server->RegisterCtor(new TIC_ctor<type>(#type,#key, getf, setf))

ChannelServer* ConstructChannelServer()
{
	ChannelServer* server = new ChannelServer();

	DEF_RW_CHANNEL(Solenoid, Bz, &Solenoid::GetBz, &Solenoid::SetBz);

	DEF_RW_CHANNEL(Quadrupole, B1, &Quadrupole::GetFieldStrength, &Quadrupole::SetFieldStrength);
	DEF_RW_CHANNEL(Sextupole, B2, &Sextupole::GetFieldStrength, &Sextupole::SetFieldStrength);
	DEF_RW_CHANNEL(Octupole, B3, &Octupole::GetFieldStrength, &Octupole::SetFieldStrength);

	DEF_RW_CHANNEL(SkewQuadrupole, B1, &SkewQuadrupole::GetFieldStrength, &SkewQuadrupole::SetFieldStrength);
	DEF_RW_CHANNEL(SkewSextupole, B2, &SkewSextupole::GetFieldStrength, &SkewSextupole::SetFieldStrength);

	DEF_RW_CHANNEL(XCor, B0, &XCor::GetFieldStrength, &XCor::SetFieldStrength);
	DEF_RW_CHANNEL(YCor, B0, &YCor::GetFieldStrength, &YCor::SetFieldStrength);

	DEF_RW_CHANNEL(SectorBend, B0, &SectorBend::GetB0, &SectorBend::SetB0);
	DEF_RW_CHANNEL(SectorBend, B1, &SectorBend::GetB1, &SectorBend::SetB1);

	DEF_RW_CHANNEL(MagnetMover, X, &MagnetMover::GetX, &MagnetMover::SetX);
	DEF_RW_CHANNEL(MagnetMover, Y, &MagnetMover::GetY, &MagnetMover::SetY);
	DEF_RW_CHANNEL(MagnetMover, Roll, &MagnetMover::GetRoll, &MagnetMover::SetRoll);

	DEF_RW_CHANNEL(TWRFStructure, E, &TWRFStructure::GetAmplitude, &TWRFStructure::SetAmplitude);
	DEF_RW_CHANNEL(TWRFStructure, Phi, &TWRFStructure::GetPhase, &TWRFStructure::SetPhase);
	DEF_RW_CHANNEL(TWRFStructure, Lambda, &TWRFStructure::GetWavelength, &TWRFStructure::SetWavelength);
	DEF_RW_CHANNEL(TWRFStructure, K, &TWRFStructure::GetK, &TWRFStructure::SetK);

	DEF_RW_CHANNEL(SWRFStructure, E, &SWRFStructure::GetAmplitude, &SWRFStructure::SetAmplitude);
	DEF_RW_CHANNEL(SWRFStructure, Phi, &SWRFStructure::GetPhase, &SWRFStructure::SetPhase);
	DEF_RW_CHANNEL(SWRFStructure, Lambda, &SWRFStructure::GetWavelength, &SWRFStructure::SetWavelength);
	DEF_RW_CHANNEL(SWRFStructure, K, &SWRFStructure::GetK, &SWRFStructure::SetK);

	DEF_RW_CHANNEL(TransverseRFStructure, E, &TransverseRFStructure::GetAmplitude,
		&TransverseRFStructure::SetAmplitude);
	DEF_RW_CHANNEL(TransverseRFStructure, Phi, &TransverseRFStructure::GetPhase, &TransverseRFStructure::SetPhase);
	DEF_RW_CHANNEL(TransverseRFStructure, Lambda, &TransverseRFStructure::GetWavelength,
		&TransverseRFStructure::SetWavelength);
	DEF_RW_CHANNEL(TransverseRFStructure, K, &TransverseRFStructure::GetK, &TransverseRFStructure::SetK);
	DEF_RW_CHANNEL(TransverseRFStructure, ROLL, &TransverseRFStructure::GetFieldOrientation,
		&TransverseRFStructure::SetFieldOrientation);

	DEF_RW_CHANNEL(CorrectorWinding, X, &CorrectorWinding::GetBy, &CorrectorWinding::SetBy);
	DEF_RW_CHANNEL(CorrectorWinding, Y, &CorrectorWinding::GetBx, &CorrectorWinding::SetBx);

	// Special BPM channels
	server->RegisterCtor(new BPMChannelCtor('X'));
	server->RegisterCtor(new BPMChannelCtor('Y'));

	return server;
}
