#include "jacobiTraits.hpp"
#include "jacobian.hpp"

int main(int argc, char **argv){

    int required=MPI_THREAD_MULTIPLE, provided;
    MPI_Init_thread(&argc, &argv, required, &provided);
    
    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    JacobianSolver L(mpi_comm);

    L.init(argc, argv);
	L.solve();

    MPI_Finalize();
};