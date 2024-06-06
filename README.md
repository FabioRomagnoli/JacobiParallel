The program can be executed as follows:
mpiexec -np $number_of_cores ./main -f $parameter_file -n $number_grid_points -t $OpenMP_threads -o $output_file

or by running the run_main.sh script that automates the command above with tweakable parameters.

The main initializes MPI and JacobiSolver.

All the parameters are either taken from command line or from param file, and recovered with getpot. Hence the class has 2 GetPot members to simplify retrieval of parameters troughout the class without storing them as members.
The following parametes can be tweaked from command line:
- number of cores (handled by MPI) 
- number of grid points
- number of OpenMP threads
- output file to save results as csv


the following from the param file:
- domain
- boundary (supports robin type of boundary. alpha and beta can be changed, and g can be a function of x and y)
- the function f
- the exact solution 
- max iterations allowed
- tolerance for convergence
- various printing and output flags further explained in the param file

If the program is called with only 1 processor, the linear solver is called, which can have OpenMP parallelization if the number of threads is > 1.

Otherwise the parallelSolver is called, all the parameters are shared between ranks, the partion is handled by the MatrixPartitioner class of the PACS examples, and from it the vectors for scatterv and gatherv as well as partition size is retrived. 
For the adjecent rows, an extra 2 rows are given to each local matrix to 
simulate the previous and sucessive block confining row, and they are sent and recieved accordingly. 
The next iteration of the solver is futher parallelized by column with the use of hybrid parallelizaion with OpenMP. Boundary conditions are then computed trough the discretization of the derivative of the neuman condition with finite differences, and appiled to the boundaries of the new iteration. 
If the solution has converged, or if the max number of iterations is reached, the loop is stopped, convergence criterias are gathered from each rank. 

The post-processing involves printing of various information if the options are activated (including exporting to vtk, or to csv). 

There is also a script run_test.sh inside the folder test, that explores the performance of the iterator with conbination of different number of cores, threads, and grid spacing. The combination wanted for the test can be changed inside the script, and it can then be executed as normal. it outputs the data in the data subfolder, and the generated plot made with a python matplotlib in the plots subfolder. 
