#include "jacobiTraits.hpp"
#include "jacobi.hpp"
#include "writeVTK.hpp"

int main(int argc, char **argv){

    int required=MPI_THREAD_MULTIPLE, provided;
    MPI_Init_thread(&argc, &argv, required, &provided);
    
    // MPI_Init(&argc, &argv);

    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    int mpi_size, mpi_rank;
    MPI_Comm_size(mpi_comm, &mpi_size);
    MPI_Comm_rank(mpi_comm, &mpi_rank);

	paramPack p;
    gridPack g;
    outputPack o; 

	Matrix f;
	Matrix exct;
	Matrix U;
    Solution sol;

    std::string outputFile; 

	if(mpi_rank == 0){
        //GetPot initialization, input takes from file param in folder
        GetPot cl(argc, argv);      //command line
        const std::string filename = cl.follow("param","-f");
        GetPot df(filename.c_str());   //datafile

        outputFile = cl.follow("output","-o");

        configOutput(cl, df, o);
	    configParams(cl, df, p, o);
        configGrid(cl, df, g, o);
        configMatrices(cl, df, g, o, f, exct, U);
	}

	if(mpi_size == 1){
        sol = linearSolver(o, p, g, f, U);
	} else {
        sol = paraSolver(mpi_comm, o,  p, g, f, U);
	}


    if(mpi_rank == 0){
	    auto [prnt_param,prnt_grid,prnt_matrix,prnt_info,prnt_result,vtk_out,csv_out] = o;
        auto [threads, maxIter,e] = p;
        auto [n,h,omega] = g;

        double err = std::sqrt(h * ((sol.Uf - exct).array().square()).sum());

        if(prnt_matrix){
            std::cout << "Solution of Jacobi Iteration" << std::endl
					 << sol.Uf << std::endl << std::endl;
        }

        if(prnt_result){
            std::cout << "The error of in the L2 norm is " << err << std::endl;
        }

        if(vtk_out){
            generateVTKFile("mesh/out.vtk", sol.Uf, n,n, h, h);
        }

        if(csv_out){
            output_dat("test/data/" + outputFile, {mpi_size, threads, n, sol.time, err});
            std::cout << "Output appended to " << outputFile + ".csv" << "\n";
        }
    }

    MPI_Finalize();

};