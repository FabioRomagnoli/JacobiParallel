#include <jacobi.hpp>

void configOutput(GetPot &cl, GetPot &df, outputPack &o){
	// Print parameters
	bool prnt_param = df("prnt_param", false);

	// Print parameters
	bool prnt_grid = df("prnt_grid", false);

	// Print matrices
	bool prnt_matrix = df("prnt_matrix", false);

	// Print info during computation
	bool prnt_info = df("prnt_info", false);

	// Print result of computation (error)
	bool prnt_result = df("prnt_result", false);

	// output to vtk for visualization
	bool vtk_out = df("vtk_out", false);

	// Output to csv for testing
	bool csv_out = df("csv_out", false);

	o = std::make_tuple(prnt_param,prnt_grid,prnt_matrix,prnt_info,prnt_result,vtk_out,csv_out);
}

void configParams(GetPot &cl, GetPot &df, paramPack &p, outputPack &o){
	auto [prnt_param,prnt_grid,prnt_matrix,prnt_info,prnt_result,vtk_out,csv_out] = o;

	// Threads
	const int threads = cl.follow(1,"-t");

	// Max iterations
    const unsigned int maxIter = df("maxIter", 1000);

	// Tollerance for convergence
    const double e = df("e", 1e-3); 

	if(prnt_param){
		std::cout << "Max iters " << maxIter << std::endl
		          << "Tolerance " << e << std::endl
				  << "Threads " << threads << std::endl << std::endl;
	}

    p = std::make_tuple(threads,maxIter,e);
}

void configGrid(GetPot &cl, GetPot &df, gridPack &g, outputPack &o){
	auto [prnt_param,prnt_grid,prnt_matrix,prnt_info,prnt_result,vtk_out,csv_out] = o;

	// Domain 
	double start = df("omega", 0.0, 0); // i specifies position if space separated array
	double end = df("omega", 1.0, 1); // i specifies position if space separated array
    const std::pair<double,double> omega = std::make_pair(start,end);

	// Number of grid points
	const int n = cl.follow(10,"-n");

	// Grid spacing
    const double h =  (omega.second - omega.first)/(n - 1);

	if(prnt_grid){
		std::cout << "Omega (" << omega.first << "," << omega.second << ")" << std::endl
				  << "Grid Points " << n << std::endl
				  << "Spacing " << h  << std::endl << std::endl;
	}

	g = std::make_tuple(n, h, omega);
}

void configMatrices(GetPot &cl, GetPot &df,  gridPack &g, outputPack &o, Matrix &f, Matrix &sol, Matrix &U){
	auto [prnt_param,prnt_grid,prnt_matrix,prnt_info,prnt_result,vtk_out,csv_out] = o;

	auto [n,h,omega] = g;

	// Boundary 
	const std::string boundary = df("boundary", "0");
	//string is parsed and a MuparserFun object is constructed from it
  	MuParserFun Fb(boundary);

	// Function 
	const std::string function = df("function", "8*(_pi^2)*sin(2*_pi*x)*sin(2*_pi*y)");
	//string is parsed and a MuparserFun object is constructed from it
  	MuParserFun Ff(function);

	// Solution 
	const std::string solution = df("solution", "sin(2*_pi*x)*sin(2*_pi*y)");
	//string is parsed and a MuparserFun object is constructed from it
  	MuParserFun Fsol(solution);

	f = Matrix::Zero(n,n); 
	sol = Matrix::Zero(n,n); 
	U = Matrix::Zero(n, n);

	for(int i = 0; i < n ; ++i){
		for(int j = 0; j < n ; ++j){
			// Imposing boundary condition
			if(i == 0 or i == n -1 or j == 0 or j == n - 1){
				f(i,j) = Fb(omega.first + i * h, omega.second + j * h);
				U(i,j) = Fb(omega.first + i * h, omega.second + j * h);
			} else {
				f(i,j) = Ff(omega.first + i * h, omega.second + j * h);
				sol(i,j) = Fsol(omega.first + i * h, omega.second + j * h);
			}
		}
	}

	if(prnt_param){
		std::cout << "Boundary " << boundary << std::endl
		          << "Function " << function << std::endl
				  << "Solution " << solution << std::endl << std::endl;
	}

	if(prnt_matrix){
		std::cout << "Function evaulated at all grid Points" << std::endl
					 << f << std::endl << std::endl 

				  << "Starting Matrix with boundary condition" << std::endl
					 << U << std::endl << std::endl 

				  << "Exact solution evaulated at all grid points" << std::endl
					 << sol << std::endl << std::endl << std::endl; 

	}
	
}