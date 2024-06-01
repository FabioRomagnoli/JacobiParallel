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

typedef std::tuple<int, unsigned int, double> paramPack;
typedef std::tuple<int, double, std::pair<double,double>> gridPack;
typedef std::tuple<bool, bool, bool,bool,  bool, bool, bool> outputPack;

struct Solution{
    Matrix Uf; 
    double time;
} ;


struct DataPoint {
    int n_cores;
    int threads;
    int n_grid_points;
    double time_elapsed;
    double error;
};

#endif /* JACOBI_TRAITS_HPP_ */