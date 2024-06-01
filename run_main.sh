#!/bin/bash

# Default values
filename="param"
num_grid_points=10
num_cores=1
threads=1

# Parse command line arguments
while getopts ":f:n:p:t:" opt; do
  case $opt in
    f)
      filename="$OPTARG"
      ;;
    n)
      num_grid_points="$OPTARG"
      ;;
    p)
      execution_mode="parallel"
      num_cores="$OPTARG"
      ;;
    t)
      num_threads="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

echo "Running main with $num_cores cores, $num_grid_points grid points, and $num_threads threads"

mpiexec -np "$num_cores" ./main -p -f "$filename" -n "$num_grid_points" -t "$num_threads"

