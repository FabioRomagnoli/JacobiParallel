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
        Solution linearSolver();

        void scatter(Matrix &lf, Matrix &lU1);
        void castAdjacentRow(Matrix &lU);
        void boundary(Matrix &U);

        Solution parallelSolver();

        void output_csv(Solution sol);


    private:
        MPI_Comm &mpi_comm;
        int mpi_size, mpi_rank;

        GetPot cl;  //command line
        GetPot df;  //datafile

        int threads;//threads
        int maxIter;
        double  tol;  //tollerance

        double alpha;
        double beta;
        Matrix g;   //boundary condition function

        int n;      //grid nodes
        double h;   //grid spacing

        apsc::MatrixPartitioner<apsc::DistributedPartitioner> mPart;
        int lrows;  //number of local rows

        Matrix f;   //function evaulated at grid
        Matrix Ue;  //exact solution
        Matrix U1;  //initial matrix

        Solution sol; 
};

#endif /* JACOBIAN_H */ 