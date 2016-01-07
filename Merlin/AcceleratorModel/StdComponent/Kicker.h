#ifndef Kicker_h
#define Kicker_h 1

#include "merlin_config.h"
// TemplateComponents
#include "AcceleratorModel/StdComponent/TemplateComponents.h"
// SimpleDrift
#include "AcceleratorModel/StdComponent/SimpleDrift.h"
// RectangularGeometry
#include "AcceleratorModel/StdGeometry/RectangularGeometry.h"

class ComponentTracker;

//	A simple thin kick section.
class Kicker : public SimpleDrift
{
public:

	explicit Kicker (const string& id, double len = 0, double x_kick = 0.0, double y_kick = 0.0, double tilt = 0.0);

	//	Returns the type string for this component.
	virtual const string& GetType () const;

	//	Returns the unique index for this class of accelerator
	//	components.
	virtual int GetIndex () const;

	//	Primary tracking interface. Prepares the specified
	//	Tracker object for tracking this component.
	virtual void PrepareTracker (ComponentTracker& aTracker);

	//	Rotates the component 180 degrees about its local Y axis.
	virtual void RotateY180 ();

	//	Virtual constructor.
	virtual ModelElement* Copy () const;

	// Data Members for Class Attributes

	//	Unique index for an Accelerator component.
	static const int ID;
	double x_kick;
	double y_kick;
	double tilt;
};


#endif
