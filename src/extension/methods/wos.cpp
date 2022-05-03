#include "wos.h"

using namespace std;

uniform_real_distribution<double> uniformDis(0.0, 1.0);
normal_distribution<double> normalDis(0,1);
const double INV_2PI = 1 / (2 * M_PI);
const double INV_4PI = (1 / (4 * M_PI));

double WOS2d::eval(double x, double y, double time) const {
	arrayd<2> p = arrayd<2>{{x,y}};
	return WOS<2>::eval(p,time);
}

double WOS2d::calcSphereArea(double r) const {
	return M_PI * r * r;
}

double WOS2d::G(arrayd<2> &x, arrayd<2> &y, double sphereR) const {
	double r = (y[0] - x[0]) * (y[0] - x[0]) + (y[1] - x[1]) * (y[1] - x[1]);
	r = sqrt(r);
	return -INV_2PI * log(sphereR / r);
}
arrayd<2> WOS2d::gradG(arrayd<2> &x, arrayd<2> &y, double sphereR) const {
	double r2 = (y[0] - x[0]) * (y[0] - x[0]) + (y[1] - x[1]) * (y[1] - x[1]);
	return (x-y)*INV_2PI*(1/r2-1/(sphereR*sphereR));
}

arrayd<2> WOS2d::randPointOnSphere(arrayd<2> &x, double radius) const {
	double theta = uniformDis(generator) * 2 * M_PI;
	return arrayd<2>{ {x[0] + radius * cos(theta), x[1] + radius * sin(theta)}};
}
arrayd<2> WOS2d::randPointOnSphere(arrayd<2> &x, double radius, arrayd<2> &opp) const {
	double theta = uniformDis(generator) * 2 * M_PI;
	opp[0] = x[0] + radius * cos(theta+M_PI);
	opp[1] = x[1] + radius * sin(theta+M_PI);
	return arrayd<2>{ {x[0] + radius * cos(theta), x[1] + radius * sin(theta)}};
}
arrayd<2> WOS2d::randPointInSphere(arrayd<2> &x, double radius) const {
	double r = radius * sqrt(uniformDis(generator));
	r = max(r, boundary->maxError);
	double theta = uniformDis(generator) * 2 * M_PI;
	return arrayd<2>{ {x[0] + r * cos(theta), x[1] + r * sin(theta)}};
}
arrayd<2> WOS2d::randPointInSphere(arrayd<2> &x, double radius, arrayd<2> &opp) const {
	double r = radius * sqrt(uniformDis(generator));
	r = max(r, boundary->maxError);
	double theta = uniformDis(generator) * 2 * M_PI;
	opp[0] = x[0] + r * cos(theta+M_PI);
	opp[1] = x[1] + r * sin(theta+M_PI);
	return arrayd<2>{ {x[0] + r * cos(theta), x[1] + r * sin(theta)}};
}



double WOS3d::calcSphereArea(double r) const {
	return M_PI * r * r * r*4/3.0;
}

double WOS3d::G(arrayd<3> &x, arrayd<3> &y, double sphereR) const {
	double r = (y-x).length();
	return -INV_4PI * (sphereR - r) / (r * sphereR);
}
arrayd<3> WOS3d::gradG(arrayd<3> &x, arrayd<3> &y, double sphereR) const {
	double r = (y-x).length();
	return (x-y)*INV_4PI*(1/(r*r*r)-1/(sphereR*sphereR*sphereR));
}

arrayd<3> WOS3d::randPointOnSphere(arrayd<3> &x, double radius) const {
	double x1 = normalDis(generator);
	double x2 = normalDis(generator);
	double x3 = normalDis(generator);
	double dist = sqrt(x1*x1 + x2*x2 + x3*x3);
	x1=(x1*radius)/dist; x2=(x2*radius)/dist; x3=(x3*radius)/dist;
	return arrayd<3> {{x.x() + x1,x.y() + x2,x.z() + x3}};
}
arrayd<3> WOS3d::randPointOnSphere(arrayd<3> &x, double radius, arrayd<3> &opp) const {
	double x1 = normalDis(generator);
	double x2 = normalDis(generator);
	double x3 = normalDis(generator);
	double dist = sqrt(x1*x1 + x2*x2 + x3*x3);
	x1=(x1*radius)/dist; x2=(x2*radius)/dist; x3=(x3*radius)/dist;
	opp[0] = x.x() - x1; opp[1] = x.y() - x2; opp[2] = x.z() - x3;
	return arrayd<3> {{x.x() + x1,x.y() + x2,x.z() + x3}};
}
arrayd<3> WOS3d::randPointInSphere(arrayd<3> &x, double radius) const {
	double ru3 = radius*cbrt(uniformDis(generator));
	ru3 = max(ru3, boundary->maxError);
	double x1 = normalDis(generator);
	double x2 = normalDis(generator);
	double x3 = normalDis(generator);

	double dist = sqrt(x1*x1 + x2*x2 + x3*x3);
	x1= (ru3*x1)/dist; x2=(ru3*x2)/dist; x3=(ru3*x3)/dist;
	return arrayd<3> {{x.x() + x1,x.y() + x2,x.z() + x3}};
}
arrayd<3> WOS3d::randPointInSphere(arrayd<3> &x, double radius, arrayd<3> &opp) const {
	double ru3 = radius*cbrt(uniformDis(generator));
	double x1 = normalDis(generator);
	double x2 = normalDis(generator);
	double x3 = normalDis(generator);

	double dist = sqrt(x1*x1 + x2*x2 + x3*x3);
	x1= (ru3*x1)/dist; x2=(ru3*x2)/dist; x3=(ru3*x3)/dist;
	opp[0] = x.x() - x1; opp[1] = x.y() - x2; opp[2] = x.z() - x3;
	return arrayd<3> {{x.x() + x1,x.y() + x2,x.z() + x3}};
}