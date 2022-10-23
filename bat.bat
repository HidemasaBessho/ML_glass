#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -M bessho@r.phys.nagoya-u.ac.jp
#$ -m ea
#$ -V
#
#$ -q all.q@einstein

./a.out $1
