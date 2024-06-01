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

	Matrix f;
	Matrix exct;
	Matrix U;
    Solution sol;

    bool vtk_out;
    bool csv_out; 

	if(mpi_rank == 0){
        //GetPot initialization, input takes from file param in folder
        GetPot command_line(argc, argv);
        const std::string filename = command_line.follow("param","-f");
        GetPot datafile(filename.c_str());

	    vtk_out = datafile("vtk_out", false);
	    csv_out = datafile("csv_out", false);

	    configParams(command_line, p, 0);
        configGrid(command_line, g);
        configMatrices(command_line,g,f,exct,U);

	}

	if(mpi_size == 1){
        sol = linearSolver(p, g, f, U);
	} else {
        sol = paraSolver(mpi_comm, p, g, f, U);
	}

    if(mpi_rank == 0){
        auto [threads, maxIter,e] = p;
        auto [n,h,omega] = g;

        double err = errL2(h,sol.Uf,exct);
        // std::cout << "Error " << err << std::endl; 

        if(vtk_out){
            generateVTKFile("mesh/out.vtk", sol.Uf, n,n, h, h);
        }

        if(csv_out){
            output_dat("test/data/output", {mpi_size, threads, n, sol.time, err});
        }
    }

    MPI_Finalize();


};