#include "jacobiOptions.hpp"
#include "jacobiTraits.hpp"

#include <iostream>

Matrix nextJacobian(const Param &p, Matrix &U){
	Matrix Up1 = Matrix::Zero(p.n,p.n);

	for(int i = 1; i < p.n - 1; ++i){
		for(int j = 1; j < p.n - 1; ++j){
			Up1(i,j) = (1.0/4.0)*(U(i-1,j) + U(i+1,j) + U(i,j-1) + 
							U(i,j+1) + (std::pow(p.h,2))*p.fEval(i,j));
		}
	}

	return Up1;
}

Matrix linearSolver(const Param &p){
	Matrix U;
	Matrix Up1 = Matrix::Zero(p.n,p.n);
	unsigned int iter = 0;
	do{		
		U = Up1;
		Up1 = nextJacobian(p, U);
		++iter;
	} while (p.e < ((Up1 - U)*p.h).norm() and iter < p.maxIter);

    return Up1;
}