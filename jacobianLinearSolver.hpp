#ifndef JACOBIANLINEARSOLVER_HPP_
#define JACOBIANLINEARSOLVER_HPP_

#include "jacobi_options.hpp"
#include "jacobi_traits.hpp"

#include <iostream>

class JacobianLinearSolver{
public:
    JacobianLinearSolver(const Param &p_):p(p_){};

    Matrix nextJacobian();
    Matrix linearSolver();

private:
    const Param &p; 
    Matrix U;
};


#endif /* JACOBIANLINEARSOLVER_HPP_ */