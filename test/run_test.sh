#!/bin/bash


paramFile="paramTest"
outputFile="test3"
cores=(1 2)
grid_points=(4 8 16)
threads=(1 2 4)

for core in "${cores[@]}"; do
  for grid in "${grid_points[@]}"; do
    for thread in "${threads[@]}"; do
      echo "Running test with $core cores, $grid grid points, and $thread threads"
      mpiexec -np "$core" ../main -f "$paramFile" -n "$grid" -t "$thread" -o "$outputFile"
    done
  done
done

python plot.py $outputFile