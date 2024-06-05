#!/bin/bash

paramFile="param"
outputFile="run"
core=3
grid=9
thread=1

mpiexec -np "$core" ./main -f "$paramFile" -n "$grid" -t "$thread" -o "$outputFile"
