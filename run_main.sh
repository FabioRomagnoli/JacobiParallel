#!/bin/bash

paramFile="param"
outputFile="run"
core=2
grid=10
thread=1

mpiexec -np "$core" ./main -f "$paramFile" -n "$grid" -t "$thread" -o "$outputFile"
