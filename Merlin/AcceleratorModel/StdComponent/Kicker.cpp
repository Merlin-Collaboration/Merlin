// Kicker
#include "AcceleratorModel/StdComponent/Kicker.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"

const int Kicker::ID = UniqueIndex();

Kicker::Kicker (const string& id, double len, double x_kick, double y_kick, double tilt)
        : x_kick(x_kick),y_kick(y_kick),tilt(tilt),TAccCompG<RectangularGeometry>(id,new RectangularGeometry(len))
{}

const string& Kicker::GetType () const
{
    _TYPESTR(Kicker);
}

int Kicker::GetIndex () const
{
    return  ID;
}

void Kicker::PrepareTracker (ComponentTracker& aTracker)
{
    _PREPTRACK(aTracker,AcceleratorComponent);
}

void Kicker::RotateY180 ()
{
    // nothing to do
}

ModelElement* Kicker::Copy () const
{
    return new Kicker(*this);
}

