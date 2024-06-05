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
        // Max iterations
        maxIter = df("maxIter", 1000);
        // Tollerance for convergence
        e = df("e", 1e-3); 

        // Domain 
        double start = df("omega", 0.0, 0); // i specifies position if space separated array
        double end = df("omega", 1.0, 1); // i specifies position if space separated array
        const std::pair<double,double> omega = std::make_pair(start,end);
        // Number of grid points
        n = cl.follow(10,"-n");
        // Grid spacing
        h =  (omega.second - omega.first)/(n - 1);

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

        U1 = Matrix::Zero(n, n);
        f = Matrix::Zero(n,n); 
        Ue = Matrix::Zero(n,n); 

        for(int i = 0; i < n ; ++i){
            for(int j = 0; j < n ; ++j){
                // Imposing boundary condition
                if(i == 0 or i == n -1 or j == 0 or j == n - 1){
                    f(i,j) = Fb(omega.first + i * h, omega.second + j * h);
                    U1(i,j) = Fb(omega.first + i * h, omega.second + j * h);
                } else {
                    f(i,j) = Ff(omega.first + i * h, omega.second + j * h);
                }
                Ue(i,j) = Fsol(omega.first + i * h, omega.second + j * h);
            }
        }

        if(df("prnt_param", false)){
            std::cout << "Max iters " << maxIter << std::endl
                    << "Tolerance " << e << std::endl
                    << "Threads " << threads << std::endl << std::endl;
            std::cout << "Boundary " << boundary << std::endl
                    << "Function " << function << std::endl
                    << "Solution " << solution << std::endl << std::endl;
            std::cout << "Omega (" << omega.first << "," << omega.second << ")" << std::endl
                    << "Grid Points " << n << std::endl
                    << "Spacing " << h  << std::endl << std::endl;
        }

        if(df("prnt_matrix", false)){
            std::cout << "Function evaulated at all grid Points" << std::endl
                        << f << std::endl << std::endl 
                    << "Starting Matrix with boundary condition" << std::endl
                        << U1 << std::endl << std::endl 
                    << "Exact solution evaulated at all grid points" << std::endl
                        << Ue << std::endl << std::endl << std::endl; 
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, mpi_comm);
    MPI_Bcast(&h, 1, MPI_DOUBLE, 0, mpi_comm);
    MPI_Bcast(&threads, 1, MPI_INT, 0, mpi_comm);
    MPI_Bcast(&maxIter, 1, MPI_INT, 0, mpi_comm);
    MPI_Bcast(&e, 1, MPI_DOUBLE, 0, mpi_comm);
}  

void JacobianSolver::solve(){
    if(mpi_size == 1){
        sol = linearSolver();
	} else {
        sol = parallelSolver();
	}

    if(mpi_rank == 0){
        if(df("prnt_matrix", false)){
            std::cout << "Solution of Jacobi Iteration" << std::endl
					 << sol.Uf << std::endl << std::endl;
        }

        if(df("prnt_result", false)){
            std::cout << "The iterator took: " << sol.time << "[s]" << std::endl;
            std::cout << "The error of in the L2 norm is " << sol.err << std::endl;
        }

        if(df("vtk_out", false)){
            generateVTKFile("mesh/out.vtk", sol.Uf, n,n, h, h);
        }

        if(df("csv_out", false)){
            output_csv(sol);
        }
    }

}

void JacobianSolver::scatter(Matrix &lf, Matrix &lU1){
    auto [counts, disp] = apsc::counts_and_displacements(mPart);

    if(mpi_rank != 0){
        f = Matrix(n,n);
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
}

void JacobianSolver::castAdjacentRow(Matrix &lU){
    /*
    Operations to send and recieve the rows necessary between adjecent proccessors to 
    compute next jacobi iteration
    */
    //each block apart from 0 sends the first row  to the block before
    if(mpi_rank > 0){
        MPI_Send(lU.row(1).data(), n, MPI_DOUBLE, mpi_rank-1, 0, mpi_comm);
    }

    //each block apart from last sends the last row to the block after
    if(mpi_rank < mpi_size - 1){
        MPI_Send(lU.row(lrows).data(), n, MPI_DOUBLE, mpi_rank+1, 1, mpi_comm);
    }
    

    //each block apart from last recieves the row from the block after
    if(mpi_rank < mpi_size - 1){
        MPI_Recv(lU.row(lrows + 1).data(), n, MPI_DOUBLE, mpi_rank+1, 0, mpi_comm, MPI_STATUS_IGNORE);
    }

    //each block apart from 0 recieves the row from the block before
    if(mpi_rank > 0){
        MPI_Recv(lU.row(0).data(), n, MPI_DOUBLE, mpi_rank-1, 1, mpi_comm, MPI_STATUS_IGNORE);
    }

    
}

void JacobianSolver::boundary(Matrix &lf, Matrix &lUk){
    //reinstate boundary conditions for the edge blocks 

    //bundary on left and right edge
#pragma omp parallel for shared(lUk, lf) num_threads(threads)
    for(int i = 0; i < lrows; ++i){
        lUk(i, 0) = lf(i, 0);
        lUk(i, n - 1) = lf(i, n - 1);
    }

    //boundary on top and bottom edge
#pragma omp parallel for shared(lUk, lf) num_threads(threads)
    for(int j = 0; j < n; ++j){
        if(mpi_rank == 0)
                lUk(0,j) = lf(0,j);

        if(mpi_rank == mpi_size - 1)
            lUk(lrows - 1 ,j) = lf(lrows - 1 ,j);
    }
}

Solution JacobianSolver::linearSolver(){
	// Print info during computation
	bool prnt_info = df("prnt_info", false);

    if(prnt_info)
        std::cout << "Problem solved linearly " << std::endl;

	Matrix Uk = U1;
    Matrix U;
    double conv;
	unsigned int iter = 0;

	tic();
	do{		
        U = Uk;
#pragma omp parallel for shared(f, U, Uk) num_threads(threads)
		for(int i = 1; i < n - 1; ++i){
			for(int j = 1; j < n - 1; ++j){
				Uk(i,j) = (1.0/4.0)*(U(i-1,j) + U(i+1,j) + U(i,j-1) + 
								U(i,j+1) + (std::pow(h,2))*f(i,j));
			}
		}
		++iter;
        conv = ((Uk - U)*h).norm();
	} while (conv > e and iter < maxIter);
	double time = toc();

    return {Uk, std::sqrt(((Ue - Uk)*h).norm()), time};
}

Solution JacobianSolver::parallelSolver(){
	// Print info during computation
	bool prnt_info = df("prnt_info", false);
    MPI_Bcast(&prnt_info, sizeof(bool), MPI_BYTE, 0, mpi_comm);

    // Partitioner from PACS examples to get vectors to give to scatter and gather
    mPart = apsc::MatrixPartitioner(n,n,mpi_size);
    auto [counts, disp] = apsc::counts_and_displacements(mPart);
    lrows = mPart.last_row(mpi_rank) - mPart.first_row(mpi_rank);

    Matrix lf(lrows, n);                     
    Matrix lU = Matrix::Zero(lrows + 2, n);  

    // Scatter to the local matrices
    scatter(lf, lU);
    
    if(prnt_info){
        if(mpi_rank == 0)
            std::cout << "Number of processes: " << mpi_size << std::endl;
        std::cout << "Number of rows on rank " << mpi_rank << ": " << lrows << std::endl;
    }

    Matrix lUk = Matrix::Zero(lrows, n);   

    double conv;
    unsigned int iter = 0; 

    tic();

    do{
#pragma omp parallel for shared(lU, lUk) num_threads(threads)
        for(int i = 0; i < lrows; ++i){
            for(int j = 0; j < n; ++j){
                lU(i + 1,j) = lUk(i,j);
            }
        }

        castAdjacentRow(lU);
 
        //compute next iteration of jacobi
#pragma omp parallel for shared(lf, lU, lUk) num_threads(threads)
        for(int i = 0; i < lrows; ++i){
            for(int j = 1; j < n  - 1; ++j){
                lUk(i,j) = 0.25*(lU(i,j) + lU(i+2,j) + 
                            lU(i+1,j-1) + lU(i+1,j+1) + (std::pow(h,2))*lf(i,j));
            }
        }

        boundary(lf,lUk);

        double lconv = (lUk - lU.middleRows(1, lrows)*h).squaredNorm();
        MPI_Allreduce(&lconv, &conv, 1, MPI_DOUBLE, MPI_SUM,mpi_comm);

        ++iter;
    }while(iter < maxIter and std::sqrt(conv) > e);

    double time = toc();

    Matrix Uf = Matrix(n,n);

    MPI_Gatherv(lUk.data(), lUk.size(), MPI_DOUBLE, Uf.data(), 
                counts.data(), disp.data(), MPI_DOUBLE, 0, mpi_comm);
    

    return {Uf, ((Ue - Uf)*h).norm(), time};
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
      << sol.time << "," << sol.err << "\n";
    
    f.close();
    std::cout << "Output appended to " << "test/data/" + outputFile + ".csv" << "\n";
}