#include "../tests.h"
#include <iostream>


#include "AcceleratorModel/Apertures/SimpleApertures.h"
#include "AcceleratorModel/Apertures/RectEllipseAperture.h"

/*
 * Test that the PointInside function for various aperture shapes
 * are working correctly.
 *
 *
*/

using namespace std;

int main(int argc, char* argv[])
{
	

	
	// Test some PointInside() functions
	RectangularAperture *rect = new RectangularAperture(12, 10);
	assert(rect->PointInside(1,1));
	assert(rect->PointInside(-1,-1));
	assert(rect->PointInside(5.9,4.9));
	assert(rect->PointInside(-5.9,4.9));
	assert(!rect->PointInside(1,5.1));
	assert(!rect->PointInside(6.1,1));
	assert(!rect->PointInside(-6.1,1));
	cout << "RectangularAperture OK" << endl;
	delete rect;

	CircularAperture *circ = new CircularAperture(12);
	assert(circ->PointInside(1,1,0));
	assert(circ->PointInside(11.9,1,0));
	assert(circ->PointInside(0, -11.9,0));
	assert(circ->PointInside(1/sqrt(2)*11.9, 1/sqrt(2)*11.9,0));
	assert(circ->PointInside(-1/sqrt(2)*11.9, -1/sqrt(2)*11.9,0));
	assert(!circ->PointInside(0, -12.1,0));
	assert(!circ->PointInside(0, 12.1,0));
	assert(!circ->PointInside(1/sqrt(2)*12, 1/sqrt(2)*12.1,0));


	cout << "CircularAperture OK" << endl;
	delete circ;

	//               double rhw, double rhh, double ehh, double ehv
	//                   RectHalfWidth, RectHalfHeight, EllipseHalfHorizontal, EllipseHalfVertical
	RectEllipseAperture *recte = new RectEllipseAperture(5, 1, 1, 2);
	assert(recte->PointInside(0.5,0.5,  0));
	assert(recte->PointInside(0.9, 0.8,  0));
	assert(recte->PointInside(0.8, 0.9,  0));
	assert(recte->PointInside(0.9, 0.1,  0));
	assert(recte->PointInside(-0.9, -0.1,  0));

	assert(!recte->PointInside(0., 1.1,  0));
	assert(!recte->PointInside(0., -1.1,  0));
	assert(!recte->PointInside(0.9, 0.9,  0));
	assert(!recte->PointInside(1, 0.1,  0));
	assert(!recte->PointInside(-1, -0.1,  0));

	cout << "RectEllipseAperture OK" << endl;
	delete recte;
}

