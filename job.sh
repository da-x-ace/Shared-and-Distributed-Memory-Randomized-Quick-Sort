#!/bin/bash

#$ -V
#$ -cwd
#$ -q development
#$ -pe 1way 36
#$ -N prob
#$ -o output_prob
#$ -e error_prob
#$ -M duke.lnmiit@gmail.com
#$ -m be
#$ -l h_rt=00:50:00

export PATH=$PATH:$HOME/cilk/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/cilk/lib
set -x
ibrun tacc_affinity ./prob /work/01905/rezaul/CSE613/HW3/tests/test-02-in.txt > test-02-out.txt
