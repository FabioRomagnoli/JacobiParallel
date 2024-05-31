#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#include <mpi.h>
#include <omp.h>
#pragma GCC diagnostic pop

#include "jacobiTraits.hpp"
#include "jacobiOptions.hpp"
#include <unistd.h>  
#include "/home/fabio/Documents/PACS/pacs-examples/Examples/src/Parallel/Utilities/partitioner.hpp"

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

static double c_start, c_diff;
#define tic() c_start = MPI_Wtime();
#define toc(x){                                      \
    c_diff = MPI_Wtime() - c_start;                  \
    std::cout << x << c_diff << " [s]" << std::endl; \
}

Matrix paraSolver(std::tuple<int,double,unsigned int,double> &pack, Matrix fEval, MPI_Comm &mpi_comm){

    int mpi_rank;
    MPI_Comm_rank(mpi_comm, &mpi_rank);

    int mpi_size;
    MPI_Comm_size(mpi_comm, &mpi_size);

    if(mpi_rank == 0)
        std::cout << "Number of processes: " << mpi_size << std::endl;

    // Set to true to print matrix, vector and result.
    const bool print = false;

    auto [n,h,maxIter,e] = pack;

    int n_rows = n;
    int n_cols = n;  
    int num_threads;
    // add openMP number of threads accordingly once implementation of openMP will be done 

    // MPI_Bcast(&num_threads, 1, MPI_INT, 0, mpi_comm);

    const unsigned int count = n_rows / mpi_size;
    const int          remainder = n_rows % mpi_size;

    const unsigned int n_rows_local =
        (mpi_rank < remainder) ? (count + 1) : count;

    std::cout << "Number of rows on rank " << mpi_rank << ": " << n_rows_local
                << std::endl;


    Matrix U = Matrix::Zero(n_rows, n_cols);
    Matrix Up1 = Matrix::Zero(n_rows, n_cols);

    //adds 2 to simplify computation and exchange of adjecent rows
    Matrix U_local = Matrix::Zero(n_rows_local + 2, n_cols);
    Matrix Up1_local = Matrix::Zero(n_rows_local, n_cols);
    Matrix fEval_local(n_rows_local, n_cols);

    // Vectors to store the number of elements to send to each
    // processor and the offset index where to start reading them from.
    std::vector<int> send_counts(mpi_size);
    std::vector<int> send_start_idx(mpi_size);

    // Vectors to store the number of elements to receive from each
    // processor and the offset index where to start writing them into.
    std::vector<int> recv_counts(mpi_size);
    std::vector<int> recv_start_idx(mpi_size);

    // MPI_Barrier(mpi_comm);
    // if(mpi_rank == 0){
    //     std::cout << fEval << std::endl;

    //     int start_idx = 0;
    //     for(int i = 0; i < mpi_size; ++i){
    //         recv_counts[i] = (i < remainder) ? (count + 1) * n_cols : count * n_cols;
    //         send_counts[i] = recv_counts[i];

    //         recv_start_idx[i] = start_idx;
    //         send_start_idx[i] = start_idx ;

    //         start_idx += recv_counts[i];
    //     }
    // }

    apsc::MatrixPartitioner mpartitioner(n_rows,n_cols,mpi_size);
    auto [counts,displacements] = apsc::counts_and_displacements(mpartitioner);

    send_counts = counts;
    send_start_idx = displacements;
    recv_counts = counts;
    recv_start_idx = displacements;

    MPI_Bcast(send_counts.data(), mpi_size, MPI_INT, 0, mpi_comm);
    MPI_Bcast(send_start_idx.data(), mpi_size, MPI_INT, 0, mpi_comm);
    MPI_Bcast(recv_counts.data(), mpi_size, MPI_INT, 0, mpi_comm);
    MPI_Bcast(recv_start_idx.data(), mpi_size, MPI_INT, 0, mpi_comm);

    MPI_Barrier(mpi_comm);
    print_vector(send_counts);
    print_vector(send_start_idx);

    MPI_Barrier(mpi_comm);
 

    tic();

    MPI_Scatterv(U.data(), send_counts.data(), send_start_idx.data(),
                MPI_DOUBLE, U_local.data() + n_cols, (n_rows_local * n_cols) , MPI_DOUBLE,
                0, mpi_comm);

    MPI_Scatterv(fEval.data(), send_counts.data(), send_start_idx.data(),
                MPI_DOUBLE, fEval_local.data(), (n_rows_local * n_cols), MPI_DOUBLE,
                0, mpi_comm);

    MPI_Barrier(mpi_comm);

   for(int rank = 0; rank < mpi_size; ++rank){
        if(rank == mpi_rank){
            std::cout << "fEval_local" << std::endl;

            std::cout << fEval_local << std::endl;
        }
        MPI_Barrier(mpi_comm);
    }

    for(int rank = 0; rank < mpi_size; ++rank){
        if(rank == mpi_rank && mpi_rank == 0){
            std::cout << std::endl;
        }
        MPI_Barrier(mpi_comm);

        if(rank == mpi_rank){
            toc("Scatter: time elapsed on rank " + std::to_string(mpi_rank) +
                ": ");
        }
        MPI_Barrier(mpi_comm);
    }

    MPI_Barrier(mpi_comm);

    int conv;
    unsigned int iter = 0; 
    do{
        std::cout << "NEW LOOP START, ITER: " << iter << " RANK: " << mpi_rank <<std::endl;

        for(int i = 0; i < n_rows_local; ++i){
            for(int j = 0; j < n_cols; ++j){
                U_local(i + 1,j) = Up1_local(i,j);
            }
        }

        /*
        Operations to send and recieve the rows necessary between adjecent proccessors to 
        compute next jacobi iteration
        */
        tic();
        //each block apart from 0 sends the first row  to the block before
        if(mpi_rank > 0){
            MPI_Send(U_local.row(1).data(), n_cols, MPI_DOUBLE, mpi_rank-1, 0, mpi_comm);
        }
        MPI_Barrier(mpi_comm);

        //each bloc apart from last sends the last row to the block after
        if(mpi_rank < mpi_size - 1){
            MPI_Send(U_local.row(n_rows_local).data(), n_cols, MPI_DOUBLE, mpi_rank+1, 1, mpi_comm);
        }
        MPI_Barrier(mpi_comm);

        //each block apart from last recieves the row from the block after
        if(mpi_rank < mpi_size - 1){
            MPI_Recv(U_local.row(n_rows_local + 1).data(), n_cols, MPI_DOUBLE, mpi_rank+1, 0, mpi_comm, MPI_STATUS_IGNORE);
        }
        MPI_Barrier(mpi_comm);

        //each block apart from 0 recieves the row from the block before
        if(mpi_rank > 0){
            MPI_Recv(U_local.row(0).data(), n_cols, MPI_DOUBLE, mpi_rank-1, 1, mpi_comm, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(mpi_comm);

        for(int rank = 0; rank < mpi_size; ++rank){
            if(rank == mpi_rank && mpi_rank == 0){
                std::cout << std::endl;
            }
            MPI_Barrier(mpi_comm);

            if(rank == mpi_rank){
                toc("Send and recieve: time elapsed on rank " + std::to_string(mpi_rank) +
                    ": ");
            }
            MPI_Barrier(mpi_comm);
        }
    
        //compute next iteration of jacobi
        for(int i = 0; i < n_rows_local; ++i){
            for(int j = 1; j < n_cols - 1; ++j){
                Up1_local(i,j) = (1.0/4.0)*(U_local(i,j) + U_local(i+2,j) + 
                            U_local(i+1,j-1) + U_local(i+1,j+1) + (std::pow(h,2))*fEval_local(i,j));
            }
        }
        MPI_Barrier(mpi_comm);

        //reinstate boundary conditions for the edge blocks 
        if(mpi_rank == 0){
            Up1_local.row(0) = Matrix::Zero(1,n_cols);
        }
        if(mpi_rank == mpi_size - 1){
            Up1_local.row(n_rows_local - 1) = Matrix::Zero(1,n_cols);
        }
        MPI_Barrier(mpi_comm);

        // Calculate the difference and apply scaling factor h
        Matrix local_diff = (Up1_local - U_local.middleRows(1, n_rows_local - 1)) * h;

        // Square each element of local_diff using array operations
        Matrix local_pow = local_diff.array().square();

        // Sum all elements of local_pow
        double local_sum = local_pow.sum();
        double total_sum; 

        MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

        // Check for convergence on rank 0 and then cast it too all ranks for loop eval
        if(mpi_rank == 0){
            bool convBool= std::sqrt(total_sum) < e;
            std::cout << std::endl << std::sqrt(total_sum) << std::endl << std::endl;
            conv = convBool ? 1 : 0;
        }
        MPI_Barrier(mpi_comm);

        MPI_Bcast(&conv, 1, MPI_INT, 0, mpi_comm);

        ++iter;

        // MPI_Barrier(mpi_comm);
        // for(int rank = 0; rank < mpi_size; ++rank){
        //     if(rank == mpi_rank){
        //         std::cout << Up1_local << std::endl;
        //     }
        //     MPI_Barrier(mpi_comm);
        // }
    

    }while(iter < maxIter and conv != 1);

    MPI_Gatherv(Up1_local.data(), Up1_local.size(), MPI_DOUBLE, Up1.data(), 
                recv_counts.data(), recv_start_idx.data(), MPI_DOUBLE, 0, mpi_comm);

    return Up1;
}