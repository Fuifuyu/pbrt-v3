#include "pbrtboundary.h"
#include "pbrt.h"

using namespace MyTypes;

double PBRTBoundary::intersectRay(arrayd<3> &x,arrayd<3> &d, arrayd<3> &n) const {
	throw "not implemented";
}

double PBRTBoundary::minDistFromBoundary(arrayd<3> &x) const {
	pbrt::Point3f p{ x.x(),x.y(),x.z() };
	double ref = distance.unsigned_distance({ p.x,p.y,p.z }).distance;
	// if (boundary->minDistanceFromPoint(p) != ref) std::cout << "minDistanceFromPoint mismatch\n";
	return ref;
}

double PBRTBoundary::minDistFromBoundary(arrayd<3> &x, arrayd<3> &d) const {
	throw "not implemented";
}