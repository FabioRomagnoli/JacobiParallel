#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "jacobiTraits.hpp"
#include "partitioner.hpp"
#include "muParserFun.hpp"
#include "writeVTK.hpp"
#include "GetPot"

#include <numbers>
#include <cmath>
#include <unistd.h>  


static double c_start, c_diff;
#define tic() c_start = MPI_Wtime();
#define toc() (c_diff = MPI_Wtime() - c_start, c_diff)


class JacobianSolver{
    public:
        JacobianSolver()= default; 
        JacobianSolver(MPI_Comm &mpi_comm_);

        void init(int argc, char **argv);

        void solve();

        void cast();
        void partition();
        void scatter(Matrix &lf, Matrix &lU1);
        void boundary(Matrix &lUk);

        Solution linearSolver();
        Solution parallelSolver();

        void output_csv(Solution sol);


    private:
        MPI_Comm &mpi_comm;
        int mpi_size, mpi_rank;

        GetPot cl;  //command line
        GetPot df;  //datafile

        int threads;//threads
        int maxIter;
        double  e;  //tollerance

        int n;      //grid nodes
        double h;   //grid spacing

        // Vectors to store the number of elements to send/recieve to/from each
        // processor and the offset index where to start reading/writing them from.
        std::vector<int> counts;
        std::vector<int> disp;

        int n_local;

        Matrix Ue;  //exact solution
        Matrix f;   //function evaulated at grid
        Matrix U1;   //initial matrix

        Solution sol; 
        
};

#endif /* JACOBIAN_H */ 