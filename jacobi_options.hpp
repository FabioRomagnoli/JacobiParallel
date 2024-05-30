#ifndef JACOBI_OPTIONS_HPP_
#define JACOBI_OPTIONS_HPP_

#include <vector>
#include "jacobi_traits.hpp"

struct Param 
{
	const F f;
    const Omega omega = Omega(0,1,0,1);
    const Boundary bounds = Boundary(0,0,0,0); 
    
    const int n = 8;
    const unsigned int n_p = 1;
    const unsigned int maxIter = 1000;

    const double h;
    const double e = 1e-3; 
    const Matrix fEval = Matrix::Zero(n,n); 
    const Matrix solEval = Matrix::Zero(n,n); 
    // const Matrix gridX = Matrix::Zero(n,n); 
    // const Matrix gridY = Matrix::Zero(n,n); 


};

#endif /* JACOBI_OPTIONS_HPP_ */