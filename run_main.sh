#!/bin/bash

paramFile="param"
outputFile="run"
core=4
grid=100
thread=1

mpiexec -np "$core" ./main -f "$paramFile" -n "$grid" -t "$thread" -o "$outputFile"
