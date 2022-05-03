#pragma once

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <functional>
#include <vector>
#include "extension/mytypes.h"

using namespace std;
using namespace MyTypes;

template <size_t N> 
struct Boundary{
    Boundary(
        function<double(arrayd<N> const&)> cond,
        double &maxError
    ) : cond(cond), maxError(maxError){};

    virtual bool isMetAt(arrayd<N> &x) const {
        return minDistFromBoundary(x) <= maxError;
    }
    virtual bool isMetAt(arrayd<N> &x, double *dist) const {
        *dist = minDistFromBoundary(x);
        return *dist < maxError;
    }
    virtual bool isInside(arrayd<N> &x) const {
        return true;
    }
    virtual double intersectRay(arrayd<N> &o, arrayd<N> &d, arrayd<N> &n) const = 0;
    virtual double minDistFromBoundary(arrayd<N> &) const = 0;
    virtual double minDistFromBoundary(arrayd<N> &, arrayd<N> &d) const = 0;

    const double &maxError;
    function<double(arrayd<N> &)> cond;
};
#endif