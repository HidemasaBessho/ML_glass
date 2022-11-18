#!/bin/sh

a=(0.56 0.72 1.0)

icc -O3 bd_N4096.cpp -o bd_N4096.out

for ((i=0 ; i<3 ; i++))
do qsub bd_N4096.bat ${a[i]}
done
