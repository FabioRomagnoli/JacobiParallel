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

class Domain{
public:
    Domain(double x0_, double xn_, double y0_, double yn_)
        : x0(x0_),xn(xn_), y0(y0_), yn(yn_) {};
    
    const double X0() const {return x0;}
    const double XN() const {return xn;}
    const double Y0() const {return y0;}
    const double YN() const {return yn;}

protected:
    const double x0; 
    const double xn; 
    const double y0; 
    const double yn;
};

class Boundary : public Domain{
public:
    Boundary(double x0_, double xn_, double y0_, double yn_)
        : Domain(x0_, xn_,y0_,yn_) {};
};

#endif /* JACOBI_TRAITS_HPP_ */