#!/bin/sh

#a=(0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.0 2.0)

icc -O3 .cpp -o .out

for ((i=0 ; i< ; i++))
do qsub .bat ${a[i]}
done
