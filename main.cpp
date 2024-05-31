#include "jacobiTraits.hpp"
#include "jacobiOptions.hpp"
#include "muParserFun.hpp"
#include "jacobiSolvers.hpp"

#include "GetPot"

#include <numbers>
#include <cmath>
#include <functional>


Scalar sol(Scalar x, Scalar y){
    return sin(2*M_PI*x)*sin(2*M_PI*y);
}

const Param getParams(GetPot &datafile, int printParam = 0){
	// Domain 
    const Omega omega = Omega(0,1,0,1);

	// Boundary condition
    const Boundary bounds = Boundary(0,0,0,0); 
    
	// Function term
	const std::string force_term = datafile("function", "8*(_pi^2)*sin(2*_pi*x)*sin(2*_pi*y)");
	//string is parsed and a MuparserFun object is constructed from it
  	MuParserFun f(force_term);

	// Number of grid points
    const int n = datafile("n", 16);

	// Grid spacing
    const double h =  (omega.XN() - omega.X0())/(n - 1);

	// Max iterations
    const unsigned int maxIter = datafile("maxIter", 1000);

	// Tollerance for convergence
    const double e = datafile("e", 1e-3); 

	// Evaluate f and exact sol at grid points once and for all
	Matrix fEval = Matrix::Zero(n,n); 
	Matrix solEval = Matrix::Zero(n,n); 

	for(int i = 1; i < n - 1; ++i){
		for(int j = 1; j < n - 1; ++j){
			fEval(i,j) = f(omega.X0()+i*h, omega.Y0()+j*h);
			solEval(i,j) = sol(omega.X0()+i*h, omega.Y0()+j*h);
		}
	}

	if(printParam){
		std::cout << "Number of grid points " << n  << std::endl
				  << "Max iterations allowed " << maxIter << std::endl
				  << "Spacing " << h << std::endl
				  << "Tolerance " << e  << std::endl;
	}

	return {omega, bounds, f, n, h, maxIter, e, fEval,solEval};
}

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);

	MPI_Comm mpi_comm = MPI_COMM_WORLD;

    int mpi_rank;
    MPI_Comm_rank(mpi_comm, &mpi_rank);

    std::tuple<int,double,unsigned int,double> pack;
	Matrix fEvalTemp;
	Matrix solEval;

	Matrix UfLinear;

	if(mpi_rank == 0){
		//GetPot initialization, input takes from file param in folder
		GetPot command_line(argc, argv);
		const std::string filename = command_line.follow("param","-f");
		GetPot datafile(filename.c_str());

		const Param p = getParams(datafile, 1);
		pack = std::tie(p.n,p.h,p.maxIter,p.e);
		fEvalTemp = p.fEval;
		solEval = p.solEval;
		//initialize and call the the linear solver
		UfLinear = linearSolver(p);
	}

	MPI_Bcast(&pack, sizeof(pack), MPI_BYTE, 0, mpi_comm);
	MPI_Barrier(mpi_comm);


	Matrix fEval = Matrix::Zero(std::get<0>(pack),std::get<0>(pack));

	if(mpi_rank == 0){
		fEval = fEvalTemp;
	}

	MPI_Bcast(fEval.data(), fEval.size(), MPI_DOUBLE, 0, mpi_comm);
    MPI_Barrier(mpi_comm);

	Matrix UfParal = paraSolver(pack, fEval, mpi_comm);
    MPI_Barrier(mpi_comm);

	if(mpi_rank == 0){
		std::cout << solEval << std::endl << std::endl;
		std::cout << UfLinear << std::endl << std::endl;
		std::cout << UfParal << std::endl << std::endl;
	}

    MPI_Finalize();
};