#include "jacobi_traits.hpp"
#include "jacobi_options.hpp"
#include "jacobianLinearSolver.hpp"

#include "GetPot"

#include <numbers>
#include <cmath>
#include <functional>

Scalar fun(Scalar x, Scalar y){
    return 8*std::pow(M_PI,2)*sin(2*M_PI*x)*sin(2*M_PI*y);
}

Scalar sol(Scalar x, Scalar y){
    return sin(2*M_PI*x)*sin(2*M_PI*y);
}

const Param initialize_params(GetPot &datafile, int printParam = 0){
	//read the function
	// const std::string fs = datafile("f", "x[0]*x[1] + 4*x[0]^(4) + x[1]^(2) + 3*x[0]");
	//parse it by giving it the number of variables present.
	// functionParser f(f);
	const F f = fun;

    const unsigned int n_p = datafile("n_p", 1);

    const unsigned int maxIter = datafile("maxIter", 1000);

	//Domain 
    const Omega omega = Omega(0,1,0,1);

	//bondary condition
    const Boundary bounds = Boundary(0,0,0,0); 
    
    const int n = datafile("n", 16);

    const double h =  (omega.getXN() - omega.getX0())/(n - 1);

    const double e = datafile("e", 1e-3); 

	//evaluate f and exact sol at grid points once and for all
	Matrix fEval = Matrix::Zero(n,n); 
	Matrix solEval = Matrix::Zero(n,n); 

	for(int i = 1; i < n - 1; ++i){
		for(int j = 1; j < n - 1; ++j){
			fEval(i,j) = f(omega.getX0()+i*h, omega.getY0()+j*h);
			solEval(i,j) = sol(omega.getX0()+i*h, omega.getY0()+j*h);

		}
	}

	//debugger matrices to check mesh is properly initialized
	// Matrix gridX = Matrix::Zero(n,n); 
	// Matrix gridY = Matrix::Zero(n,n); 
	// for(int i = 0; i < n; ++i){
	// 	for(int j = 0; j < n; ++j){
	// 		gridX(i,j) = omega.getX0()+i*h;
	// 		gridY(i,j) = omega.getY0()+j*h;
	// 	}
	// }

	if(printParam){
		std::cout << "Number of processors " << n_p  << std::endl
				<< "Max iterations allowed " << maxIter << std::endl
				<< "Number of grid points " << n  << std::endl
				<< "Tolerance " << e  << std::endl
				<< "Spacing " << h << std::endl;
	}

	return {f, omega, bounds, n, n_p, maxIter, h, e, fEval, solEval};
}

int main(int argc, char **argv){
	//GetPot initialization, input takes from file param in folder
	GetPot command_line(argc, argv);
  	const std::string filename = command_line.follow("param","-f");
  	GetPot datafile(filename.c_str());

	const Param p = initialize_params(datafile, 1);

	//initialize and call the the linear solver
	JacobianLinearSolver jLinSolver(p);
	Matrix Uf = jLinSolver.linearSolver();

	std::cout << p.solEval << std::endl;
};