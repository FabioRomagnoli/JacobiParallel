#!/bin/bash

filename="param"
cores=(1 2 4)
grid_points=(4 8 16 32 64)
threads=(1 2)

for core in "${cores[@]}"; do
  for grid in "${grid_points[@]}"; do
    for thread in "${threads[@]}"; do
      echo "Running test with $core cores, $grid grid points, and $thread threads"
      mpiexec -np "$core" ./main -p -f "$filename" -n "$grid" -t "$thread"
    done
  done
done


python ./test/data/plot.py