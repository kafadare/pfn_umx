#!/bin/bash
#$ -N bivar_umx
#$ -o /cbica/projects/bgd-pfn/bivar_umx_output/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -j y

#submit this with something like the following command
#qsub -t 1:n cubicscriptdraft.sh

source /cbica/projects/bgd-pfn/.miniconda3/etc/profile.d/conda.sh

conda activate umx-4.15

home=/cbica/projects/bgd-pfn/
echo $SGE_TASK_ID
INDEX=$SGE_TASK_ID  #this is a special variable that is replaced with 1:n
echo $INDEX

Rscript  ${home}bivar_pfn.R $INDEX
