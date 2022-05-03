#pragma once

#ifndef PDES_H
#define PDES_H

#include "extension/core/pde.h"

struct TaylorVortexPDE2 : PDE<2>{
    TaylorVortexPDE2():PDE(PDETypes::Poisson){};

    virtual double laplacianOP(arrayd<2> &x, double t) override;
    virtual double truth(arrayd<2> &x, double t) override;
    virtual arrayd<2> gradTruth(arrayd<2> & x, double t) override;
};

struct TaylorVortexPDE3 : PDE<3>{
    TaylorVortexPDE3():PDE(PDETypes::Poisson){};

    virtual double laplacianOP(arrayd<3> &x, double t) override;
    virtual double truth(arrayd<3> &x, double t) override;
    virtual arrayd<3> gradTruth(arrayd<3> & x, double t) override;
};
#endif