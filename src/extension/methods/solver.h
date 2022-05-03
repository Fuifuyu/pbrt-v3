#pragma once

#ifndef SOLVER_H
#define SOLVER_H

#include "extension/mytypes.h"

template <size_t N>
struct Boundary;
template <size_t N>
struct PDE;

using namespace MyTypes;

template <size_t N>
class Solver{
    public:
        Solver(){};
        Solver(PDE<N> *pde, Boundary<N> *boundary):pde{pde},boundary{boundary}{};
        virtual arrayd<N> evalGrad(arrayd<N> &x, double t) const = 0;
        virtual arrayd<N> evalGrad(arrayd<N> &x, double t, std::vector<double> &u) const = 0;
        virtual double eval(arrayd<N> &x, double t) const = 0;
        virtual double eval(arrayd<N> &x, double t, std::vector<double> &u) const = 0;
        virtual double cond(arrayd<N> &x) const {
            return boundary->cond(x);
        }
        double truth(arrayd<N> &x, double t) const{
            return pde->truth(x,t);
        }
        arrayd<N> gradTruth(arrayd<N> &x, double t) const{
            return pde->gradTruth(x,t);
        }
    protected:
        PDE<N> *pde;
        Boundary<N> *boundary;
};

#endif