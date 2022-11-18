#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q@einstein

./bd_N4096.out $1
