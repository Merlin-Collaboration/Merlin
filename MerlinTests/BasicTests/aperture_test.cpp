#include "../tests.h"
#include <iostream>


#include "SimpleApertures.h"
#include "RectEllipseAperture.h"
#include "CollimatorAperture.h"

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
	cout << "RectangularAperture ";
	cout.flush();
	RectangularAperture *rect = new RectangularAperture(12, 10);
	assert(rect->PointInside(1,1));
	assert(rect->PointInside(-1,-1));
	assert(rect->PointInside(5.9,4.9));
	assert(rect->PointInside(-5.9,4.9));
	assert(!rect->PointInside(1,5.1));
	assert(!rect->PointInside(6.1,1));
	assert(!rect->PointInside(-6.1,1));
	cout << "OK" << endl;
	delete rect;

	cout << "CircularAperture ";
	cout.flush();
	CircularAperture *circ = new CircularAperture(12);
	assert(circ->PointInside(1,1,0));
	assert(circ->PointInside(11.9,1,0));
	assert(circ->PointInside(0, -11.9,0));
	assert(circ->PointInside(1/sqrt(2)*11.9, 1/sqrt(2)*11.9,0));
	assert(circ->PointInside(-1/sqrt(2)*11.9, -1/sqrt(2)*11.9,0));
	assert(!circ->PointInside(0, -12.1,0));
	assert(!circ->PointInside(0, 12.1,0));
	assert(!circ->PointInside(1/sqrt(2)*12, 1/sqrt(2)*12.1,0));
	cout << "OK" << endl;
	delete circ;

	//               double rhw, double rhh, double ehh, double ehv
	//                   RectHalfWidth, RectHalfHeight, EllipseHalfHorizontal, EllipseHalfVertical
	cout << "RectEllipseAperture ";
	cout.flush();
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

	cout << "OK" << endl;
	delete recte;

	/*  For checking try the following in python
	 *  t = numpy.linspace(0,2*pi, 500)
	 *  rx = numpy.array([-1,-1,1,1,-1])
	 *  ry = numpy.array([-1,1,1,-1,-1])
	 *  pyplot.plot(5*rx, 1*ry)
	 *  pyplot.plot(1*sin(t), 2*cos(t))
	 *  pyplot.grid()
	 *  pyplot.show()
	 */




	cout << "CollimatorAperture ";
	cout.flush();
	// double w,double h, double t, Material* m, double length, double x_off, double y_off
	CollimatorAperture *collapp = new CollimatorAperture(12, 10, 0, nullptr, 1, 0, 0);
	collapp->SetExitWidth(collapp->GetFullEntranceWidth());
	collapp->SetExitHeight(collapp->GetFullEntranceHeight());
	assert(collapp->PointInside(1,1,0));
	assert(collapp->PointInside(-1,-1,0));
	assert(collapp->PointInside(5.9,4.9,0));
	assert(collapp->PointInside(-5.9,4.9,0));
	assert(!collapp->PointInside(1,5.1,0));
	assert(!collapp->PointInside(6.1,1,0));
	assert(!collapp->PointInside(-6.1,1,0));

	assert(collapp->PointInside(1,1,1));
	assert(collapp->PointInside(-1,-1,1));
	assert(collapp->PointInside(5.9,4.9,1));
	assert(collapp->PointInside(-5.9,4.9,1));
	assert(!collapp->PointInside(1,5.1,1));
	assert(!collapp->PointInside(6.1,1,1));
	assert(!collapp->PointInside(-6.1,1,1));

	cout << "OK" << endl;
	delete collapp;


	cout << "CollimatorAperture with tilt ";
	collapp = new CollimatorAperture(24, 20, M_PI/4, nullptr, 1, 0, 0);
	collapp->SetExitWidth(collapp->GetFullEntranceWidth());
	collapp->SetExitHeight(collapp->GetFullEntranceHeight());
	assert(collapp->PointInside(1,1,0));
	assert(collapp->PointInside(-1,-1,0));
	assert(collapp->PointInside(1.4,15.5,0));
	assert(!collapp->PointInside(1.3,15.5,0));
	assert(!collapp->PointInside(1.5,15.5,0));
	assert(!collapp->PointInside(1.4,15.6,0));
	assert(collapp->PointInside(-10.2,-6.7,0));
	assert(!collapp->PointInside(-10.2,-6.8,0));
	assert(collapp->PointInside(15.5,1.4,0));
	assert(collapp->PointInside(14.0,-0.1,0));
	assert(!collapp->PointInside(14.0,-0.2,0));

	assert(collapp->PointInside(1,1,1));
	assert(collapp->PointInside(-1,-1,1));
	assert(collapp->PointInside(1.4,15.5,1));
	assert(!collapp->PointInside(1.3,15.5,1));
	assert(!collapp->PointInside(1.5,15.5,1));
	assert(!collapp->PointInside(1.4,15.6,1));
	assert(collapp->PointInside(-10.2,-6.7,1));
	assert(!collapp->PointInside(-10.2,-6.8,1));
	assert(collapp->PointInside(15.5,1.4,1));
	assert(collapp->PointInside(14.0,-0.1,1));
	assert(!collapp->PointInside(14.0,-0.2,1));


	cout << "OK" << endl;
	delete collapp;


	/*  For checking try the following in python
	 *  rx = numpy.array([-1,-1,1,1,-1])
	 *  ry = numpy.array([-1,1,1,-1,-1])
	 *  pyplot.plot(*zip(*[ rotate(x,y, pi/4) for x,y in zip(rx*12, ry*10) ]))
	 *  pyplot.grid()
	 *  pyplot.show()
	 */

	cout << "CollimatorAperture with offset";
	collapp = new CollimatorAperture(2, 2, 0, nullptr, 1, 1, 1);
	collapp->SetExitWidth(collapp->GetFullEntranceWidth());
	collapp->SetExitHeight(collapp->GetFullEntranceHeight());
	collapp->SetExitXOffset(collapp->GetEntranceXOffset());
	collapp->SetExitYOffset(collapp->GetEntranceYOffset());
	assert(collapp->PointInside(1,1,0));
	assert(collapp->PointInside(0.1,0.1,0));
	assert(!collapp->PointInside(-0.1,0.1,0));
	assert(!collapp->PointInside(0.1,-0.1,0));
	assert(collapp->PointInside(1.9,1.9,0));
	assert(!collapp->PointInside(2.1,1.9,0));
	assert(!collapp->PointInside(1.9,2.1,0));
	assert(!collapp->PointInside(2.1,2.1,0));

	assert(collapp->PointInside(1,1,1));
	assert(collapp->PointInside(0.1,0.1,1));
	assert(!collapp->PointInside(-0.1,0.1,1));
	assert(!collapp->PointInside(0.1,-0.1,1));
	assert(collapp->PointInside(1.9,1.9,1));
	assert(!collapp->PointInside(2.1,1.9,1));
	assert(!collapp->PointInside(1.9,2.1,1));
	assert(!collapp->PointInside(2.1,2.1,1));

	cout << "OK" << endl;
	delete collapp;


	cout << "CollimatorAperture with taper";
	collapp = new CollimatorAperture(4, 4, 0, nullptr, 1, 0, 0);
	collapp->SetExitWidth(2);
	collapp->SetExitHeight(2);
	collapp->SetExitXOffset(collapp->GetEntranceXOffset());
	collapp->SetExitYOffset(collapp->GetEntranceYOffset());
	assert(collapp->PointInside(1.9,1.9,0));
	assert(!collapp->PointInside(1.9,1.9,1));
	assert(collapp->PointInside(-1.9,1.9,0));
	assert(!collapp->PointInside(-1.9,1.9,1));
	assert(collapp->PointInside(1.9,-1.9,0));
	assert(!collapp->PointInside(1.9,-1.9,1));

	assert(collapp->PointInside(1.4,1.4,0.5));
	assert(!collapp->PointInside(1.6,1.6,0.5));
	assert(collapp->PointInside(-1.4,1.4,0.5));
	assert(!collapp->PointInside(-1.6,1.6,0.5));
	assert(collapp->PointInside(-1.4,-1.4,0.5));
	assert(!collapp->PointInside(-1.6,-1.6,0.5));

	cout << "OK" << endl;
	delete collapp;


}

