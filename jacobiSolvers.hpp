#ifndef JACOBISOLVERS_H
#define JACOBISOLVERS_H

#include "jacobiTraits.hpp"
#include "jacobiOptions.hpp"

Matrix nextJacobian(const Param &p, Matrix &U);
Matrix linearSolver(const Param &p);

Matrix paraSolver(std::tuple<int,double,unsigned int,double> &pack, Matrix fEval, MPI_Comm &mpi_comm);


#endif /* JACOBISOLVERS_H */ 
