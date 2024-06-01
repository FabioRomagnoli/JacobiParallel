#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#include <mpi.h>
#include <omp.h>
#pragma GCC diagnostic pop

#include <jacobi.hpp>

Solution paraSolver(MPI_Comm &mpi_comm, outputPack &o, paramPack &p, gridPack &g, Matrix &f, Matrix &U){

    int mpi_size, mpi_rank;
    MPI_Comm_size(mpi_comm, &mpi_size);
    MPI_Comm_rank(mpi_comm, &mpi_rank);

    MPI_Bcast(&o, sizeof(o), MPI_BYTE, 0, mpi_comm);
    MPI_Bcast(&p, sizeof(p), MPI_BYTE, 0, mpi_comm);
    MPI_Bcast(&g, sizeof(g), MPI_BYTE, 0, mpi_comm);

	auto [prnt_param,prnt_grid,prnt_matrix,prnt_info,prnt_result,vtk_out,csv_out] = o;
    auto [threads, maxIter,e] = p;
    auto [n,h,omega] = g;

    const bool print = false;


    MPI_Barrier(mpi_comm);
    if(mpi_rank == 0 and prnt_info)
        std::cout << "Number of processes: " << mpi_size << std::endl;

    // Set to true to print matrix, vector and result.

    int n_rows = n;
    int n_cols = n;  
    int num_threads = threads; 


    // Vectors to store the number of elements to send/recieve to/from each
    // processor and the offset index where to start reading/writing them from.
    std::vector<int> counts(mpi_size);
    std::vector<int> displacements(mpi_size);

    // Partitioner from PACS examples to get vectors to give to scatter and gather
    apsc::MatrixPartitioner mpartitioner(n_rows,n_cols,mpi_size);
    auto count_disp = apsc::counts_and_displacements(mpartitioner);

    counts = count_disp[0];
    displacements = count_disp[1];

    const unsigned int n_rows_local = mpartitioner.last_row(mpi_rank) - mpartitioner.first_row(mpi_rank);

    MPI_Barrier(mpi_comm);
    
    if(prnt_info)
        std::cout << "Number of rows on rank " << mpi_rank << ": " << n_rows_local << std::endl;

    //adds 2 rows to simplify computation and exchange of adjecent rows
    Matrix U_local = Matrix::Zero(n_rows_local + 2, n_cols);
    Matrix Up1_local = Matrix::Zero(n_rows_local, n_cols);
    Matrix f_local(n_rows_local, n_cols);
    Matrix Up1 = Matrix::Zero(n_rows, n_cols);

    if(mpi_rank != 0){
        f = Matrix::Zero(n,n); 
	    U = Matrix::Zero(n,n); 
    }

    MPI_Barrier(mpi_comm);
    MPI_Scatterv(U.data(), counts.data(), displacements.data(),
                MPI_DOUBLE, U_local.data() + n_cols, (n_rows_local * n_cols) , MPI_DOUBLE,
                0, mpi_comm);

    MPI_Barrier(mpi_comm);
    MPI_Scatterv(f.data(), counts.data(), displacements.data(),
                MPI_DOUBLE, f_local.data(), (n_rows_local * n_cols), MPI_DOUBLE,
                0, mpi_comm);

    int conv;
    unsigned int iter = 0; 

    tic();

    do{

#pragma omp parallel for shared(U_local, Up1_local) num_threads(num_threads)
        for(int i = 0; i < n_rows_local; ++i){
            for(int j = 0; j < n_cols; ++j){
                U_local(i + 1,j) = Up1_local(i,j);
            }
        }

        /*
        Operations to send and recieve the rows necessary between adjecent proccessors to 
        compute next jacobi iteration
        */
        //each block apart from 0 sends the first row  to the block before
        if(mpi_rank > 0){
            MPI_Send(U_local.row(1).data(), n_cols, MPI_DOUBLE, mpi_rank-1, 0, mpi_comm);
        }

        //each block apart from last sends the last row to the block after
        if(mpi_rank < mpi_size - 1){
            MPI_Send(U_local.row(n_rows_local).data(), n_cols, MPI_DOUBLE, mpi_rank+1, 1, mpi_comm);
        }
        MPI_Barrier(mpi_comm);

        //each block apart from last recieves the row from the block after
        if(mpi_rank < mpi_size - 1){
            MPI_Recv(U_local.row(n_rows_local + 1).data(), n_cols, MPI_DOUBLE, mpi_rank+1, 0, mpi_comm, MPI_STATUS_IGNORE);
        }

        //each block apart from 0 recieves the row from the block before
        if(mpi_rank > 0){
            MPI_Recv(U_local.row(0).data(), n_cols, MPI_DOUBLE, mpi_rank-1, 1, mpi_comm, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(mpi_comm);

        //compute next iteration of jacobi
#pragma omp parallel for shared(f_local, U_local, Up1_local) num_threads(num_threads)
        for(int i = 0; i < n_rows_local; ++i){
            for(int j = 1; j < n_cols - 1; ++j){
                Up1_local(i,j) = (1.0/4.0)*(U_local(i,j) + U_local(i+2,j) + 
                            U_local(i+1,j-1) + U_local(i+1,j+1) + (std::pow(h,2))*f_local(i,j));
            }
        }

        //reinstate boundary conditions for the edge blocks 
        if(mpi_rank == 0){
            for(int j = 0; j < n_cols ; ++j){
                Up1_local(0,j) = f(0,j);
            }
        }
        if(mpi_rank == mpi_size - 1){
            for(int j = 0; j < n_cols ; ++j){
                Up1_local(n_rows_local - 1,j) = f(n_rows_local - 1 ,j);
            }
        }


        // Calculate the difference and apply scaling factor h
        Matrix local_diff = (Up1_local - U_local.middleRows(1, n_rows_local)) * h;

        // Square each element of local_diff using array operations
        Matrix local_pow = local_diff.array().square();

        // Sum all elements of local_pow
        double local_sum = local_pow.sum();
        double total_sum; 

        MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

        // Check for convergence on rank 0 and then cast it too all ranks for loop eval
        if(mpi_rank == 0){
            bool convBool= std::sqrt(total_sum) < e;
            // std::cout << std::endl << std::sqrt(total_sum) << std::endl << std::endl;
            conv = convBool ? 1 : 0;
        }
        MPI_Barrier(mpi_comm);


        MPI_Bcast(&conv, 1, MPI_INT, 0, mpi_comm);
        MPI_Barrier(mpi_comm);

        ++iter;
    }while(iter < maxIter and conv != 1);

    double time_elapsed = toc();

    MPI_Gatherv(Up1_local.data(), Up1_local.size(), MPI_DOUBLE, Up1.data(), 
                counts.data(), displacements.data(), MPI_DOUBLE, 0, mpi_comm);
    
    return {Up1, time_elapsed};
}