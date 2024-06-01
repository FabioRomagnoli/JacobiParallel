#include "jacobi.hpp"

Solution linearSolver(outputPack &o, paramPack &p, gridPack &g, Matrix &f, Matrix &U){
	auto [prnt_param,prnt_grid,prnt_matrix,prnt_info,prnt_result,vtk_out,csv_out] = o;
    auto [threads, maxIter,e] = p;
    auto [n,h,omega] = g;

    if(prnt_info)
        std::cout << "Problem solved linearly " << std::endl;

	Matrix Up1 = Matrix::Zero(n,n);
	unsigned int iter = 0;

	tic();

	do{		
		U = Up1;
#pragma omp parallel for shared(f, U, Up1) num_threads(threads)
		for(int i = 1; i < n - 1; ++i){
			for(int j = 1; j < n - 1; ++j){
				Up1(i,j) = (1.0/4.0)*(U(i-1,j) + U(i+1,j) + U(i,j-1) + 
								U(i,j+1) + (std::pow(h,2))*f(i,j));
			}
		}
		++iter;
	} while (e < ((Up1 - U)*h).norm() and iter < maxIter);

	double time_elapsed = toc();

    return {Up1, time_elapsed};
}