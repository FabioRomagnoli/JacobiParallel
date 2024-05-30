#include "jacobianLinearSolver.hpp"

Matrix JacobianLinearSolver::nextJacobian(){
	Matrix Up1 = Matrix::Zero(p.n,p.n);

	for(int i = 1; i < p.n - 1; ++i){
		for(int j = 1; j < p.n - 1; ++j){
			Up1(i,j) = (1.0/4.0)*(U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) + (std::pow(p.h,2))*p.fEval(i,j));
		}
	}

	return Up1;
}

Matrix JacobianLinearSolver::linearSolver(){
	Matrix Up1 = Matrix::Zero(p.n,p.n);
	idx iter = 0;
	do{		
		U = Up1;
		Up1 = nextJacobian();
		++iter;
	} while (p.e < ((Up1 - U)*p.h).norm() and iter < p.maxIter);

	std::cout << U - p.solEval << std::endl << std::endl;

    return Up1;
}