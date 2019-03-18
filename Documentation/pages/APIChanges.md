# API Changes {#APIChanges}

# Changes by version {#changesbyversion}

This page documents changes to the to MERLIN that may require updates to
user programs.

[TOC]

## Version 5.02 {#APIChanges502}

### Behaviour change in LatticeFunctionTable

Behaviour of calling LatticeFunctionTable::Calculate() has changed in the case where the longitudinal motion is unstable (for example lattices without RF). Previously it would return a table, but some of the values would be invalid (NaN). Now it will throw a MerlinException if you try to calculate a function using non transverse terms. If you were previously checking the content of the table for NaNs, you should now use a try/catch block to catch the error. Also see LatticeFunctionTable::SetForceLongitudinalStability() which can simplify the process of dealing with lattices without RF.

### Deprecation of ScaleBendPathLength and RingDeltaTProcess

For lattices with no RF the longitudinal part of the transfer matrix will be unstable, which prevents LatticeFunctionTable calculations from working. Previously by setting a bend scale, RingDeltaTProcess was used to adjust the ct particle coordinate make the dynamics stable. This however required finding a value of the bend scale that would work for a given lattice.

This has now been replaced by a simpler option.

LatticeFunctionTable::SetForceLongitudinalStability();

This achieves the same by directly changing the R_55 transfer matrix element.

### MADInterface changes

The MADInterface class has been refactored to utilize the newly developed and versitile DataTable container class. Moreover, the ConstructApertures function has been removed as it was a redundant function designed for direct aperture construction from a MAD input file. The function was also notably incomplete and buggy. All intended functionality is still provided via the ApertureConfiguration class, which has been also been redesigned for simplicity and modularity.

### ApertureConfiguration changes

The ConfigureElementApertures function has been streamlined to utilise the new Aperture function's improved factory pattern design. As a result, any non-viable aperture input will return an error and the default aperture is automatically a RECTELLIPSE of dimensions {1,1,1,1}. This has made the SetDefaultAperture and EnableDefaultAperture member function obsolete and they have been removed, accordingly.

### RandomNG changes

RandomNG has switched from using the internal ACG random number generator to the Mersenne Twister algorithm found in the C++11 standard library. There are only minor changes to the RandomNG api.

Mersenne Twister can take a vector of unsigned ints as the seed. This allows more flexible uses of seeds without risk of collision. For example study involving multiple parameter sets and multiple cluster jobs might use a seed for each. The multiple seeds can be set for example like:

    RandomNG::init();
    RandomNG::init(seed);// where seed is a std::uint32_t
    RandomNG::init(seeds); // where seeds is a std::vector<std::uint32_t>
    RandomNG::init({seed1, seed2});

The generator will also be self seeded (using the random values from the operating system), if RandomNG::init() is not called or called with no arguments.

Likewise RandomNG::getSeed() will now return a vector:

    static std::vector<std::uint32_t> getSeed();

RandGenerator is also removed, code using it to get access to its own generator should instead use RandomNG::getLocalGenerator().

Remove all classes related to old random number generators: ACG, Binomial, DiscUnif, Erlang, Geom, HypGeom, Landau, LogNorm, MLCG, NegExp, Normal, Poisson, RNG, Random, RndInt, Uniform, Weibull.

Some additional functions have been added, see documentation for RandomNG.

### Removal of ParticleBunchConstructor

The ParticleBunchConstructor has been removed and its functionality moved into the constructor of the ParticleBunch class. This simplifies creating bunches and makes it easy to add new distribution types.

For example the previous use of ParticleBunchConstructor:

    ParticleBunchConstructor* constructor = new ParticleBunchConstructor(mybeam, npart, normalDistribution());
    ProtonBunch* myBunch = constructor->ConstructParticleBunch<ProtonBunch>();

can be replaced by:

    ProtonBunch* myBunch = new ProtonBunch(npart, NormalParticleDistributionGenerator(), mybeamm);

Where the distribution type was previously set by an enumeration, it is now set by passing a ParticleDistributionGenerator. The equivalents are:

    normalDistribution              ->   NormalParticleDistributionGenerator()
    flatDistribution                ->   UniformParticleDistributionGenerator()
    ringDistribution                ->   RingParticleDistributionGenerator()
    skewHaloDistribution            ->   RingParticleDistributionGenerator()
    horizontalHaloDistribution1     ->   HorizonalHalo1ParticleDistributionGenerator()
    verticalHaloDistribution1       ->   VerticalHalo1ParticleDistributionGenerator()
    horizontalHaloDistribution2     ->   HorizonalHalo2ParticleDistributionGenerator()
    verticalHaloDistribution2       ->   VerticalHalo2ParticleDistributionGenerator()

Instead of using ParticleBunchConstructor::SetFilter(), filters can be passed into ParticleBunch:

    ParticleBunch myBunch(npart, UniformParticleDistributionGenerator(), mybeamm, myfilter);

Additional options beyond the scope of BeamData can be passed as arguments to the ParticleDistributionGenerator. For example 3 sigma cuts in the normal distribution:

    NormalParticleDistributionGenerator(3.0);

See ParticleTracking::ParticleBunch, ParticleDistributionGenerator

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
