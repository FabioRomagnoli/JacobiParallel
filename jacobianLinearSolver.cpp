#include "jacobiTraits.hpp"
#include "jacobi.hpp"

Matrix linearSolver(int argc, char **argv){

	paramPack p;
	Matrix fEval;
	Matrix solEval;
	Matrix U;

	configParams(argc, argv, p, fEval, solEval, U, 0);
    auto [n,h,maxIter,e,omega] = p;

	Matrix Up1 = Matrix::Zero(n,n);
	unsigned int iter = 0;

	do{		
		U = Up1;
		for(int i = 1; i < n - 1; ++i){
			for(int j = 1; j < n - 1; ++j){
				Up1(i,j) = (1.0/4.0)*(U(i-1,j) + U(i+1,j) + U(i,j-1) + 
								U(i,j+1) + (std::pow(h,2))*fEval(i,j));
			}
		}
		++iter;
	} while (e < ((Up1 - U)*h).norm() and iter < maxIter);

    return Up1;
}