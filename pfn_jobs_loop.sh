#!/bin/bash
#$ -N pfn_loop
#$ -o /cbica/projects/bgd-pfn/$JOB_NAME.log
#$ -j y

for i in {1..2}
do
for k in {1..2}
export network=$i
export vertex=$k
qsub umx_on_vertex.sh -V
sleep 30
done
sleep 60
done
