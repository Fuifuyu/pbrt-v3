#define _USE_MATH_DEFINES
#include <math.h>
#include "extension/pde/pdes.h"
double const MIN_X = 3*M_PI/2;
double const PI_3 = 2*M_PI;

double TaylorVortexPDE2::laplacianOP(arrayd<2> & x, double t){
    return 8 * M_PI * M_PI * sin(PI_3*x.x()) * sin(PI_3*x.y());
}

double TaylorVortexPDE2::truth(arrayd<2> & x, double t){
    return -sin(PI_3*x.x()) * sin(PI_3*x.y());
}

arrayd<2> TaylorVortexPDE2::gradTruth(arrayd<2> & x, double t){
    return {{-PI_3*cos(PI_3*x.x())*sin(PI_3*x.y()),-PI_3*sin(PI_3*x.x())*cos(PI_3*x.y())}};
}


double TaylorVortexPDE3::laplacianOP(arrayd<3> & x, double t){
    return 12 * M_PI * M_PI * sin(PI_3*x.x()) * sin(PI_3*x.y()) * sin(PI_3*x.z());
}

double TaylorVortexPDE3::truth(arrayd<3> & x, double t){
    double s1 = sin(PI_3 * x.x());
    double s2 = sin(PI_3 * x.y());
    double s3 = sin(PI_3 * x.z());
    return -s1 * s2 * s3;
}

arrayd<3> TaylorVortexPDE3::gradTruth(arrayd<3> & x, double t){
    return {{
        -PI_3*cos(PI_3*x.x())*sin(PI_3*x.y())*sin(PI_3*x.z()),
        -PI_3*sin(PI_3*x.x())*cos(PI_3*x.y())*sin(PI_3*x.z()),
        -PI_3*sin(PI_3*x.x())*sin(PI_3*x.y())*cos(PI_3*x.z()),
        }};
}