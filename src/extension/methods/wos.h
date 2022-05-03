#pragma once
#define _USE_MATH_DEFINES

#ifndef WOS_H
#define WOS_H

#include "extension/core/pde.h"
#include "extension/mytypes.h"
#include "extension/core/boundary.h"
#include "solver.h"
#include <array>
#include <random>
#include <math.h>

using namespace std;
using namespace MyTypes;

template<size_t N>
class WOS : public Solver<N>{
public:
	WOS(){};
	WOS(PDE<N> *pde, Boundary<N> *boundary, mt19937 generator)
	: Solver{pde,boundary}, generator(generator)
	{};
	virtual double eval(arrayd<N> &x, double t) const override
	{
		double R;
		if (boundary->isMetAt(x, &R)) return boundary->cond(x);
		arrayd<N> randOnSphere = randPointOnSphere(x, R);
		arrayd<N> randInSphere = randPointInSphere(x, R);

		double u_1 = eval(randOnSphere, t);
		double harmonic = 0;
		if (pde->type == PDETypes::Poisson) {
			double f_1 = pde->laplacianOP(randInSphere, t);
			harmonic = f_1 * G(x, randInSphere, R) * calcSphereArea(R);
		}
		return u_1 + harmonic;
	}
	virtual double eval(arrayd<N> &x, double t, std::vector<double> &u) const override{
		throw exception("Not implemented");
	};
	virtual arrayd<N> evalGrad(arrayd<N> &x, double t) const override{
		double R = boundary->minDistFromBoundary(x);
		arrayd<N> opp;
		arrayd<N> randOnSphere = randPointOnSphere(x, R, opp);
		arrayd<N> sol = (randOnSphere-x)*(N*eval(randOnSphere,t)*0.5/(R*R));
		sol += (opp-x)*(N*eval(opp,t)*0.5/(R*R));
		if(pde->type==PDETypes::Poisson){
			arrayd<N> randInSphere = randPointInSphere(x,R, opp);
			arrayd<N> source = gradG(x,randInSphere,R)*calcSphereArea(R)*pde->laplacianOP(randInSphere,t)*0.5;
			 source += gradG(x,opp,R)*M_PI*R*R*pde->laplacianOP(opp,t)*0.5;
			return sol + source;
		}
		return sol;
	}
    virtual arrayd<N> evalGrad(arrayd<N> &x, double t, std::vector<double> &u) const override{
		throw exception("Not implemented");
	};
protected:
	mutable std::mt19937 generator;

	virtual double calcSphereArea(double radius) const = 0;
	virtual arrayd<N> randPointOnSphere(arrayd<N> &x, double radius) const = 0;
	virtual arrayd<N> randPointOnSphere(arrayd<N> &x, double radius, arrayd<N> &opp) const = 0;
	virtual arrayd<N> randPointInSphere(arrayd<N> &x, double radius) const = 0;
	virtual arrayd<N> randPointInSphere(arrayd<N> &x, double radius, arrayd<N> &opp) const = 0;
	virtual double G(arrayd<N> &x, arrayd<N> &y, double sphereR) const = 0;
	virtual arrayd<N> gradG(arrayd<N> &x, arrayd<N> &y, double sphereR) const = 0;
};

class WOS2d : public WOS<2> {
public:
	using WOS::WOS;
	using WOS::eval;

	double eval(double x, double y, double time) const;
protected:
	double calcSphereArea(double radius) const override;
	arrayd<2> randPointOnSphere(arrayd<2> &x, double radius) const override;
	arrayd<2> randPointInSphere(arrayd<2> &x, double radius) const override;
	arrayd<2> randPointOnSphere(arrayd<2> &x, double radius, arrayd<2> &opp) const override;
	arrayd<2> randPointInSphere(arrayd<2> &x, double radius, arrayd<2> &opp) const override;
	double G(arrayd<2> &x, arrayd<2> &y, double sphereR) const override;
	arrayd<2> gradG(arrayd<2> &x, arrayd<2> &y, double sphereR) const override;
};

class WOS3d : public WOS<3> {
public:
	using WOS::WOS;

protected:
	double calcSphereArea(double radius) const override;
	arrayd<3> randPointOnSphere(arrayd<3> &x, double radius) const override;
	arrayd<3> randPointInSphere(arrayd<3> &x, double radius) const override;
	arrayd<3> randPointOnSphere(arrayd<3> &x, double radius, arrayd<3> &opp) const override;
	arrayd<3> randPointInSphere(arrayd<3> &x, double radius, arrayd<3> &opp) const override;
	double G(arrayd<3> &x, arrayd<3> &y, double sphereR) const override;
	arrayd<3> gradG(arrayd<3> &x, arrayd<3> &y, double sphereR) const override;
};

#endif