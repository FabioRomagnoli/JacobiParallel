#include "jacobian.hpp"

JacobianSolver::JacobianSolver(MPI_Comm &mpi_comm_):mpi_comm(mpi_comm_){
    MPI_Comm_size(mpi_comm, &mpi_size);
    MPI_Comm_rank(mpi_comm, &mpi_rank);
}  

void JacobianSolver::init(int argc, char **argv){
    if(mpi_rank == 0){
        //GetPot initialization, input takes from file param in folder
        cl = GetPot(argc, argv);      //command line
        const std::string filename = cl.follow("param","-f");
        df = GetPot(filename.c_str());   //datafile

        // OpenMP threads
        threads = cl.follow(1,"-t");

        // Domain 
        double start = df("domain/omega", 0.0, 0); // i specifies position if space separated array
        double end = df("domain/omega", 1.0, 1); // i specifies position if space separated array
        const std::pair<double,double> omega = std::make_pair(start,end);

        // Number of grid points
        n = cl.follow(10,"-n");
        // Grid spacing
        h =  (omega.second - omega.first)/(n - 1);

        // Max iterations
        maxIter = df("parameters/maxIter", 10000);
        // Tolerance for convergence
        tol = df("parameters/tol", 1e-3); 

        // Boundary coefficients
        alpha = df("boundary/alpha", 1.0);
        beta = df("boundary/beta", 0.5);

        // Boundary 
        const std::string Sg = df("boundary/g", "0");
        MuParserFun Fg(Sg);

        // Function 
        const std::string function = df("problem/function", "8*(_pi^2)*sin(2*_pi*x)*sin(2*_pi*y)");
        MuParserFun Ff(function);

        // Solution 
        const std::string solution = df("problem/solution", "0");
        //string is parsed and a MuparserFun object is constructed from it
        MuParserFun Fsol(solution);

        f = Matrix::Zero(n,n); 
        g = Matrix::Zero(n,n); 
        U1 = Matrix::Zero(n, n);
        Ue = Matrix(n,n); 

        for(int i = 0; i < n ; ++i){
            for(int j = 0; j < n ; ++j){
                Ue(i,j) = Fsol(omega.first + i * h, omega.first + j * h);

                if(i == 0 or i == n -1 or j == 0 or j == n - 1){
                    g(i,j) = Fg(omega.first + i * h, omega.first + j * h);
                } else {
                    f(i,j) = Ff(omega.first + i * h, omega.first + j * h);
                }
            }
        }

        boundary(U1);

        if(df("print/param", false)){
            std::cout 
                    << "    Info:  "  << std::endl
                    << "Cores = " << mpi_size << std::endl
                    << "Threads = " << threads << std::endl << std::endl

                    << "Domain = (" << omega.first << "," << omega.second << ")" << std::endl
                    << "Grid Points = " << n << std::endl
                    << "Spacing = " << h  << std::endl << std::endl 

                    << "Max iters = " << maxIter << std::endl
                    << "Tolerance = " << tol << std::endl << std::endl

                    << "    Boundary:   " << std::endl
                    << "alpha = " << alpha << std::endl
                    << "beta = " << beta << std::endl
                    << "Boundary: g(x,y) = " << Sg << std::endl << std::endl

                    << "    Functions:   " << std::endl
                    << "Function: f(x,y) = " << function << std::endl
                    << "Solution: u(x,y) = " << solution << std::endl << std::endl;
        }

        if(df("print/matrix", false)){
            std::cout 
                    << "Function evaluated at all grid points" << std::endl
                        << f << std::endl << std::endl;
                
            if(!Ue.isZero()){
                std::cout << "Exact solution evaluated at all grid points" << std::endl
                        << Ue << std::endl << std::endl; 
            } else {
                std::cout << "The exact solution is not known" << std::endl << std::endl;  
            }
                   
        }
    }

    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, mpi_comm);
    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, mpi_comm);
    MPI_Bcast(&n, 1, MPI_INT, 0, mpi_comm);
    MPI_Bcast(&h, 1, MPI_DOUBLE, 0, mpi_comm);
    MPI_Bcast(&threads, 1, MPI_INT, 0, mpi_comm);
    MPI_Bcast(&maxIter, 1, MPI_INT, 0, mpi_comm);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, mpi_comm);
}  

void JacobianSolver::solve(){
    if(mpi_size == 1){
        sol = linearSolver();
	} else {
        sol = parallelSolver();
	}

    if(mpi_rank == 0){
        if(df("print/matrix", false)){
            std::cout << "Solution with Jacobi Iteration" << std::endl
					 << sol.Uf << std::endl << std::endl;
        }

        if(df("print/result", false)){
            std::cout << "The iterator took: " << sol.time << "[s]" << std::endl;
            std::cout << "The iterations converged to: c = " << sol.conv << std::endl;
            if(!Ue.isZero()){
                std::cout << "The error in the L2 norm: e = " << ((Ue - sol.Uf)*h).norm() << std::endl;
            } else {
                std::cout << "The exact solution is not known" << std::endl;
            }
        }

        if(df("output/vtk_out", false)){
            generateVTKFile("mesh/out.vtk", sol.Uf, n,n, h, h);
        }

        if(df("output/csv_out", false)){
            output_csv(sol);
        }
    }

}

Solution JacobianSolver::linearSolver(){
	// Print info during computation
	bool info = df("print/info", false);

    if(info)
        std::cout << "Problem solved linearly " << std::endl;

	Matrix Uk(n,n);
    Matrix U = U1;
    double conv;
	unsigned int iter = 0;

	tic();
	do{		
#pragma omp parallel for shared(f, U, Uk) num_threads(threads) collapse(2)
		for(int i = 1; i < n - 1; ++i){
			for(int j = 1; j < n - 1; ++j){
				Uk(i,j) = (1.0/4.0)*(U(i-1,j) + U(i+1,j) + U(i,j-1) + 
								U(i,j+1) + (std::pow(h,2))*f(i,j));
			}
		}
		++iter;
        conv = ((Uk - U)*h).norm();
        U = Uk;
	} while (conv > tol and iter < maxIter);
	double time = toc();
    
    return {U, conv, time};
}

void JacobianSolver::scatter(Matrix &lf, Matrix &lU1){
    auto [counts, disp] = apsc::counts_and_displacements(mPart);

    if(mpi_rank != 0){
        f = Matrix(n,n);
        g = Matrix(n,n);
	    U1 = Matrix(n,n);
        Ue = Matrix(n,n);
    }

    MPI_Scatterv(U1.data(), counts.data(), disp.data(),
                MPI_DOUBLE, lU1.data() + n, (lrows * n) , MPI_DOUBLE,
                0, mpi_comm);

    MPI_Scatterv(f.data(), counts.data(), disp.data(),
                MPI_DOUBLE, lf.data(), (lrows * n), MPI_DOUBLE,
                0, mpi_comm);

    MPI_Bcast(Ue.data(), Ue.size(), MPI_DOUBLE, 0, mpi_comm);
    MPI_Bcast(g.data(), g.size(), MPI_DOUBLE, 0, mpi_comm);

}

void JacobianSolver::castAdjacentRow(Matrix &lU){
    /* Operations to send and recieve the rows necessary between adjecent proccessors to 
    compute next jacobi iteration */

    //each block apart from 0 sends the first row  to the block before
    if(mpi_rank > 0)
        MPI_Send(lU.row(1).data(), n, MPI_DOUBLE, mpi_rank-1, 0, mpi_comm);

    //each block apart from last sends the last row to the block after
    if(mpi_rank < mpi_size - 1)
        MPI_Send(lU.row(lrows).data(), n, MPI_DOUBLE, mpi_rank+1, 1, mpi_comm);

    //each block apart from last recieves the row from the block after
    if(mpi_rank < mpi_size - 1)
        MPI_Recv(lU.row(lrows + 1).data(), n, MPI_DOUBLE, mpi_rank+1, 0, mpi_comm, MPI_STATUS_IGNORE);
    
    //each block apart from 0 recieves the row from the block before
    if(mpi_rank > 0)
        MPI_Recv(lU.row(0).data(), n, MPI_DOUBLE, mpi_rank-1, 1, mpi_comm, MPI_STATUS_IGNORE);
}

void JacobianSolver::boundary(Matrix &U) {
    // Update boundary points using the Robin boundary condition
    int rows = U.rows();
    int cols = U.cols();

    double denom = alpha * h + beta;

    for (int i = 0; i < rows; ++i) {
        // Left boundary (j = 0)
        U(i, 0) = (h * g(i, 0) + beta * U(i, 1)) / (denom);

        // Right boundary (j = n-1)
        U(i, cols - 1) = (h * g(i, cols - 1) + beta * U(i, cols - 2)) / (denom);
    }
    
    for(int j = 0; j < cols; ++j){
        if(mpi_rank == 0)
            U(0,j) = (h * g(0,j) + beta * U(1,j)) / (denom);

        if(rows < n){
            // Top boundary (i = n - 1)
            if(mpi_rank == mpi_size - 1)
                U(rows - 1,j) = (h * g(n - 1,j) + beta * U(rows - 2,j)) / (denom);
        } else {
            U(rows - 1,j) = (h * g(rows - 1,j) + beta * U(rows - 2,j)) / (denom);
        }
        
    }
}

Solution JacobianSolver::parallelSolver(){
	// Print info during computation
	bool info = df("print/info", false);
    MPI_Bcast(&info, sizeof(bool), MPI_BYTE, 0, mpi_comm);

    // Partitioner from PACS examples to get vectors to give to scatter and gather
    mPart = apsc::MatrixPartitioner(n,n,mpi_size);
    auto [counts, disp] = apsc::counts_and_displacements(mPart);
    lrows = mPart.last_row(mpi_rank) - mPart.first_row(mpi_rank);

    Matrix lf(lrows, n);                     
    Matrix lU = Matrix::Zero(lrows + 2, n);  

    // Scatter to the local matrices
    scatter(lf, lU);
    
    if(info){
        std::cout << "Number of rows on rank " << mpi_rank << ": " << lrows << std::endl;
    }

    Matrix lUk = Matrix::Zero(lrows, n);   

    double conv;
    unsigned int iter = 0; 

    tic();

    do{
        castAdjacentRow(lU);
 
        //compute next iteration of jacobi
#pragma omp parallel for shared(lf, lU, lUk) num_threads(threads) collapse(2)
        for(int i = 0; i < lrows; ++i){
            for(int j = 1; j < n - 1; ++j){
                lUk(i,j) = 0.25*(lU(i,j) + lU(i+2,j) + 
                            lU(i+1,j-1) + lU(i+1,j+1) + (std::pow(h,2))*lf(i,j));
            }
        }

        boundary(lUk);

        double lconv = ((lUk - lU.middleRows(1, lrows))*h).squaredNorm();
        MPI_Allreduce(&lconv, &conv, 1, MPI_DOUBLE, MPI_SUM,mpi_comm);

#pragma omp parallel for shared(lU, lUk) num_threads(threads) collapse(2)
        for(int i = 0; i < lrows; ++i){
            for(int j = 0; j < n; ++j){
                lU(i + 1,j) = lUk(i,j);
            }
        }

        ++iter;
    }while(iter < maxIter and std::sqrt(conv) > tol);

    double time = toc();

    Matrix Uf = Matrix(n,n);

    MPI_Gatherv(lUk.data(), lUk.size(), MPI_DOUBLE, Uf.data(), 
                counts.data(), disp.data(), MPI_DOUBLE, 0, mpi_comm);
    

    return {Uf, std::sqrt(conv), time};
}

void JacobianSolver::output_csv(Solution sol) {
    std::string outputFile = cl.follow("output","-o");
    std::ofstream f;
    f.open("data/" + outputFile + ".csv", std::ios::app);

    if (!f.is_open()) {
        std::cerr << "Error: Could not open file " << outputFile + ".csv" << std::endl;
        return;
    }

    // Check if file is empty to write header
    f.seekp(0, std::ios::end);
    if (f.tellp() == 0) {
        // Write the header
        f << "n_cores,threads,n_grid_points,time_elapsed,error\n";
    }

    f << mpi_size << "," << threads << "," << n << ","
      << sol.time << "," << ((Ue - sol.Uf)*h).norm() << "\n";
    
    f.close();
    std::cout << "Output appended to " << "test/data/" + outputFile + ".csv" << "\n";
}
