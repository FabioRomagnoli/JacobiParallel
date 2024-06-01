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


void configParams(GetPot &datafile, paramPack &p, int printParam);

void configGrid(GetPot &datafile, gridPack &g);

void configMatrices(GetPot &datafile, gridPack &g, Matrix &f, Matrix &sol, Matrix &U);

void print_vector(const std::vector<int>& vec);

double errL2(double h, Matrix &Uf, Matrix &sol);

Solution linearSolver(paramPack &p, gridPack &g, Matrix &f, Matrix &U);

Solution paraSolver(MPI_Comm &mpi_comm, paramPack &p, gridPack &g, Matrix &f, Matrix &U);

void output_dat(const std::string& filename, DataPoint dp);


#endif /* JACOBI_H */ 
