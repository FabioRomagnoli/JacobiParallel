#include <jacobi.hpp>

void configParams(int argc, char **argv, paramPack &p, Matrix &fEval, Matrix &solEval, Matrix &U, int printParam = 1){
	//GetPot initialization, input takes from file param in folder
    GetPot command_line(argc, argv);
    const std::string filename = command_line.follow("param","-f");
    GetPot datafile(filename.c_str());

	// Domain 
    const std::pair<double,double> omega = std::make_pair(0.0,1.0);

	// Number of grid points
    const int n = datafile("n", 16);

	// Grid spacing
    const double h =  (omega.second - omega.first)/(n - 1);

	// Max iterations
    const unsigned int maxIter = datafile("maxIter", 1000);

	// Tollerance for convergence
    const double e = datafile("e", 1e-3); 

	if(printParam){
		std::cout << "Number of grid points " << n  << std::endl
				  << "Max iterations allowed " << maxIter << std::endl
				  << "Spacing " << h << std::endl
				  << "Tolerance " << e  << std::endl;
	}

    p = std::make_tuple(n,h,maxIter,e,omega);

	// Boundary condition
    const Boundary bounds = Boundary(0,0,0,0);

	// Function 
	const std::string function = datafile("function", "8*(_pi^2)*sin(2*_pi*x)*sin(2*_pi*y)");
	//string is parsed and a MuparserFun object is constructed from it
  	MuParserFun f(function);

	// Solution 
	const std::string solution = datafile("solution", "sin(2*_pi*x)*sin(2*_pi*y)");
	//string is parsed and a MuparserFun object is constructed from it
  	MuParserFun sol(solution);

	fEval = Matrix::Zero(n,n); 
	solEval = Matrix::Zero(n,n); 

	for(int i = 1; i < n - 1; ++i){
		for(int j = 1; j < n - 1; ++j){
			fEval(i,j) = f(omega.first + i * h, omega.second + j * h);
			solEval(i,j) = sol(omega.first + i * h, omega.second + j * h);
		}
	}

	U = Matrix::Zero(n, n);
}




void print_vector(const std::vector<int>& vec) {
    std::cout << std::endl;
    std::cout << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]";
    std::cout << std::endl;
}