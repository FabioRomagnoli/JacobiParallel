#ifndef JACOBI_TRAITS_HPP_
#define JACOBI_TRAITS_HPP_
#include "Eigen/Core"

#include <functional>
#include <math.h>

using Scalar = double;
using idx = size_t;
using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
using F = std::function<Scalar(Scalar, Scalar)>;


class Domain{
public:
    Domain(double x0_, double xn_, double y0_, double yn_)
        : x0(x0_),xn(xn_), y0(y0_), yn(yn_) {};
    
    const double getX0() const {return x0;}
    const double getXN() const {return xn;}
    const double getY0() const {return y0;}
    const double getYN() const {return yn;}

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

class Omega : public Domain{
public:
    Omega(double x0_, double xn_, double y0_, double yn_)
        : Domain(x0_, xn_,y0_,yn_) {};
};

#endif /* JACOBI_TRAITS_HPP_ */