#!/bin/bash


filename="paramTest"
output="test2"
cores=(1 2 4)
grid_points=(4 8 16 32 64 128 256)
threads=(1 2)

for core in "${cores[@]}"; do
  for grid in "${grid_points[@]}"; do
    for thread in "${threads[@]}"; do
      echo "Running test with $core cores, $grid grid points, and $thread threads"
      mpiexec -np "$core" ./main -f "$filename" -n "$grid" -t "$thread" -o "$output"
    done
  done
done

python test/plot.py $output