#include "jacobiTraits.hpp"
#include "jacobi.hpp"

int main(int argc, char **argv){
	Matrix Usol;

	if(true){
		Usol = paraSolver(argc, argv);
	} else {
		Usol = linearSolver(argc, argv);
	}
};