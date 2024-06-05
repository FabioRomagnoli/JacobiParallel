#ifndef JACOBI_TRAITS_HPP_
#define JACOBI_TRAITS_HPP_
#include "Eigen/Core"
#include <Eigen/Dense>
#include "mpi.h"
#include <math.h>


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <tuple>


using Scalar = double;
using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using F = std::function<Scalar(Scalar, Scalar)>;

struct Solution{
    Matrix Uf; 
    double err;
    double time;
} ;


#endif /* JACOBI_TRAITS_HPP_ */