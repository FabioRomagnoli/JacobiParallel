#ifndef JACOBI_H
#define JACOBI_H

#include "jacobiTraits.hpp"
#include "partitioner.hpp"
#include "muParserFun.hpp"
#include "GetPot"

#include <numbers>
#include <cmath>
#include <unistd.h>  


static double c_start, c_diff;
#define tic() c_start = MPI_Wtime();
#define toc() (c_diff = MPI_Wtime() - c_start, c_diff)


// Functions to configure the iterator
void configOutput(GetPot &cl, GetPot &df, outputPack &o);

void configParams(GetPot &cl, GetPot &df, paramPack &p, outputPack &o);

void configGrid(GetPot &cl, GetPot &df, gridPack &g, outputPack &o);

void configMatrices(GetPot &cl, GetPot &df, gridPack &g, outputPack &o, Matrix &f, Matrix &sol, Matrix &U);


// Functions to solve jacobi iteration linearly or parallelized 
Solution linearSolver(outputPack &o, paramPack &p, gridPack &g, Matrix &f, Matrix &U);

Solution paraSolver(MPI_Comm &mpi_comm, outputPack &o,  paramPack &p, gridPack &g, Matrix &f, Matrix &U);


// Utility functions to print vectors and output data to file
void output_dat(const std::string& filename, DataPoint dp);

void print_vector(const std::vector<int>& vec);


#endif /* JACOBI_H */ 
