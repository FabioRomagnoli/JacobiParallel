#include <jacobi.hpp>

void configParams(GetPot &command_line, paramPack &p, int printParam = 0){
	const std::string filename = command_line.follow("param","-f");
	GetPot datafile(filename.c_str());

	// Threads
	const int threads = command_line.follow(1,"-t");


	// Max iterations
    const unsigned int maxIter = datafile("maxIter", 1000);

	// Tollerance for convergence
    const double e = datafile("e", 1e-3); 

	if(printParam){
		std::cout << "Max iterations allowed " << maxIter << std::endl
				  << "Tolerance " << e  << std::endl;
	}

    p = std::make_tuple(threads,maxIter,e);
}

void configGrid(GetPot &command_line, gridPack &g){
	// Domain 
    const std::pair<double,double> omega = std::make_pair(0.0,1.0);

	// Number of grid points
	const int n = command_line.follow(10,"-n");

	// Grid spacing
    const double h =  (omega.second - omega.first)/(n - 1);

	g = std::make_tuple(n, h, omega);
}

void configMatrices(GetPot &command_line,  gridPack &g, Matrix &f, Matrix &sol, Matrix &U){
	const std::string filename = command_line.follow("param","-f");
	GetPot datafile(filename.c_str());

	auto [n,h,omega] = g;

	// Boundary condition
    const Boundary bounds = Boundary(0,0,0,0);

	// Function 
	const std::string function = datafile("function", "8*(_pi^2)*sin(2*_pi*x)*sin(2*_pi*y)");
	//string is parsed and a MuparserFun object is constructed from it
  	MuParserFun Ff(function);

	// Solution 
	const std::string solution = datafile("solution", "sin(2*_pi*x)*sin(2*_pi*y)");
	//string is parsed and a MuparserFun object is constructed from it
  	MuParserFun Fsol(solution);

	f = Matrix::Zero(n,n); 
	sol = Matrix::Zero(n,n); 

	for(int i = 1; i < n - 1; ++i){
		for(int j = 1; j < n - 1; ++j){
			f(i,j) = Ff(omega.first + i * h, omega.second + j * h);
			sol(i,j) = Fsol(omega.first + i * h, omega.second + j * h);
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

double errL2(double h, Matrix &Uf, Matrix &sol){
	return std::sqrt(h * ((Uf - sol).array().square()).sum());
}

void output_dat(const std::string& filename, DataPoint dp) {
    std::ofstream f;
    f.open(filename + ".csv", std::ios::app);

    if (!f.is_open()) {
        std::cerr << "Error: Could not open file " << filename + ".csv" << std::endl;
        return;
    }

    // Check if file is empty to write header
    f.seekp(0, std::ios::end);
    if (f.tellp() == 0) {
        // Write the header
        f << "n_cores,threads,n_grid_points,time_elapsed,error\n";
    }

    f << dp.n_cores << "," << dp.threads << "," << dp.n_grid_points << ","
      << dp.time_elapsed << "," << dp.error << "\n";
    
    f.close();
    std::cout << "Output appended to " << filename + ".csv" << "\n";
}