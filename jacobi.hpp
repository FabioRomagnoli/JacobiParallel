#ifndef JACOBI_H
#define JACOBI_H

#include "jacobiTraits.hpp"

#include "muParserFun.hpp"
#include "GetPot"
#include <numbers>
#include <cmath>
#include <functional>
#include <unistd.h>  
#include "partitioner.hpp"


void configParams(int argc, char **argv, paramPack &p, Matrix &fEval, Matrix &solEval, Matrix &U, int printParam);

void print_vector(const std::vector<int>& vec);

Matrix linearSolver(int argc, char **argv);

Matrix paraSolver(int argc, char **argv);

#endif /* JACOBI_H */ 
