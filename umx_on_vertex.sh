#!/bin/sh
#$ -N pfn_vertex_umx
#$ -o /cbica/projects/bgd-pfn/pfn_vertex_output/$JOB_NAME_$JOB_ID.log
#$ -j y
# Set the amount of memory being requested.
#$ -l s_vmem=10G

#submit this with something like the following command
#qsub -t 1:n cubicscriptdraft.sh
source /cbica/projects/bgd-pfn/.miniconda3/etc/profile.d/conda.sh
conda activate umx-4.15

home=/cbica/projects/bgd-pfn/pfn_umx/
cd $home
#network=$i
#vertex=$k
#vertex=$SGE_TASK_ID  #this is a special variable that is replaced with 1:n

#matlab -r "call_matlab_fct($INDEX, '/cbica/projects/bgd-pfn/pfn1_vertex_pairs.csv');quit"
matlab -r "get_vertex_col($network,$vertex);quit"
Rscript  ${home}vertex_umx.R $network $vertex
