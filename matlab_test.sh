#!/bin/sh
#$ -N matlab
#$ -o /cbica/projects/bgd-pfn/pfn_vertex_output/$JOB_NAME_$JOB_ID.log
#$ -j y
# Set the amount of memory being requested.
#$ -l h_vmem=8G

# Source the project-user’s miniconda setup
source /cbica/projects/bgd-pfn/.miniconda3/etc/profile.d/conda.sh
conda activate umx-4.15
# Input arguments
n=$1
k=$2
#OUTDIR=$3

# Load and unload modules needed to run your code
#module load matlab
# Run matlab script
#matlab -r "call_matlab_fct($1, '/cbica/projects/bgd-pfn/pfn1_vertex_pairs.csv');quit"
matlab -r "get_vertex_col($n,$k);quit"
