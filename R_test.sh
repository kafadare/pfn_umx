#!/bin/sh
#
# The name of the job
#$ -N matlab_test
#
#
# Set the amount of memory being requested.
#$ -l h_vmem=8G

# Source the project-user’s miniconda setup

#activate conda base environments
conda activate base

# Input arguments
n=$1
k=$2
#OUTDIR=$3

# Load and unload modules needed to run your code
#module load matlab
# Run R script
Rscript vertex_umx.R '/cbica/projects/bgd-pfn/vertex_columns/V$n_PFN$k.csv' > vertex_umx_test.txt
