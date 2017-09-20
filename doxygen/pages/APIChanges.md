# API Changes {#APIChanges}

# Changes by version {#changesbyversion}

This page documents changes to the to MERLIN that may require updates to
user programs.

[TOC]

## Version 5.02 {#APIChanges502}

### Indecies -> Indexes

The spelling of several function names was changed:

    AcceleratorComponent::AppendBeamlineIndecies() ->  AcceleratorComponent::AppendBeamlineIndexes()
    AcceleratorModel::GetIndecies()                ->  AcceleratorModel::GetIndexes()
    CollimateParticleProcess::GetIndecies()        ->  CollimateParticleProcess::GetIndexes()
    ComponentFrame::AppendBeamlineIndecies()       ->  ComponentFrame::AppendBeamlineIndexes()
    CorrectorWinding::AppendBeamlineIndecies()     ->  CorrectorWinding::AppendBeamlineIndexes()
    Klystron::AppendBeamlineIndecies()             ->  Klystron::AppendBeamlineIndexes()
    ModelElement::GetBeamlineIndecies()            ->  ModelElement::GetBeamlineIndexes()
    SequenceFrame::AppendBeamlineIndecies()        ->  SequenceFrame::AppendBeamlineIndexes()

See AcceleratorComponent, AcceleratorModel, ParticleTracking::CollimateParticleProcess, 
ComponentFrame, CorrectorWinding, Klystron, ModelElement, SequenceFrame

## Version 5.01 {#APIChanges501}

### Directory Flattening

Subdirectories within the source code have been removed, so includes need
to be updated. e.g.

    #include "NumericalUtils/PhysicalUnits.h"

to

    #include "PhysicalUnits.h"

This can be done automatically with the command:

    sed -i 's@#include "[A-Za-z/]*/@#include "@g' filename.cpp

### Dustbin -> CollimationOutput

The Dustbin classes have been renamed to ParticleTracking::CollimationOutput .

    LossMapDustbin* myLossMapDustbin = new LossMapDustbin(tencm);
    myCollimateProcess->SetDustbin(myLossMapDustbin);

to

    LossMapCollimationOutput* myLossOutput = new LossMapCollimationOutput(tencm);
    myCollimateProcess->SetCollimationOutput(myLossOutput);

Likewise the header file is renamed:

    Dustbin.h -> LossMapCollimationOutput.h

See ParticleTracking::CollimationOutput.

### ScatteringModel

Rather than using Collimation::ScatteringModel and calling SetScatterType(), you can
now use a specific preconfigured model, e.g. ScatteringModelMerlin

See Collimation::ScatteringModel, Collimation::ScatteringModelMerlin, Collimation::ScatteringModelSixTrack, Collimation::ScatteringModelSixTrackElastic, Collimation::ScatteringModelSixTrackIoniz, Collimation::ScatteringModelSixTrackSD.
