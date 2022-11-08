#!/bin/sh
#$ -N pfn_loop
#$ -o /cbica/projects/bgd-pfn/$JOB_NAME.log
#$ -j y

for i in {1..2}
do
for k in {1..2}
qsub umx_on_vertex.sh -F "network=$i, vertex=$k"
sleep 30
done
sleep 60
done
