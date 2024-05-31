#ifndef JACOBI_OPTIONS_HPP_
#define JACOBI_OPTIONS_HPP_

#include "jacobiTraits.hpp"

struct Param 
{    
    Omega omega = Omega(0,1,0,1);
    Boundary bounds = Boundary(0,0,0,0); 
    F f;

    int n = 8;
    double h;

    unsigned int maxIter = 1000;
    double e = 1e-3; 
    
    Matrix fEval = Matrix::Zero(n,n); 
    Matrix solEval = Matrix::Zero(n,n); 
};

#endif /* JACOBI_OPTIONS_HPP_ */