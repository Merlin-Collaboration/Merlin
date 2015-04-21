#include "AcceleratorModel/Aperture.h"

ostream& operator<< (ostream& out, const Aperture& ap) {
    ap.printout(out);
    return out;
}

void Aperture::printout(std::ostream& out) const{
	out << GetApertureType();
}

