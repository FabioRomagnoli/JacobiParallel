The program can be executed as follows:
mpiexec -np $number_of_cores ./main -f $parameter_file -n $number_grid_points -t $OpenMP_threads -o $output_file

The main initializes a few variables of type tuple defined in jacobiTraits to simplify the process of communicating them to the various processors. They take information both from command line: 
- number of cores (handled by MPI) 
- number of grid points
- number of OpenMP threads
- output file to save results
and from the param file:
- domain
- boundary (controleld by a function of x and y and computed only at the boundary)
- the function to f
- the exact solution 
- max iterations allowed
- tolerance for convergence
- various printing flags further explained in the param file

They are all initialized with their respecive config function only in rank 0, to be then casted to the ranks if and as needed. 

If the program is called with only 1 processor, the linear solver is called, which can have OpenMP parallelization if the number of threads is > 1.

Otherwise the parallelSolver is called, all the parameters are shared between ranks, the partion is handled by the MatrixPartitioner class of the PACS examples, and from it the vectors for scatterv and gatherv as well as partition size is retrived. For the adjecent rows, an extra 2 rows are given to each local matrix to 
simulate the previous and sucessive block confining row, and they are sent and recieved accordingly. The next iteration of the solver is futher parallelized by column with the use of hybrid parallelizaion with OpenMP. Boundary conditions are reinstated for edge blocks. If the matrix has converged, or if the max number of iterations is reached, the loop is stopped and results are gathered. 

The post-processing involves printing of various information if the options are activated (including exporting to vtk, or to csv). 

There is also a script run_test.sh, that explores the performance of the iterator with different number of cores, threads, and grid spacing. it just has to be executed as a normal script (it contains various options to refine the testing) and outputs the data in the data folder and the generated comparison plot made with a python script using matplotlib in the plots folder. 
