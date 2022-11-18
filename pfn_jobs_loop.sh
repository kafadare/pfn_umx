#!/bin/bash
#$ -N pfn_loop
#$ -o /cbica/projects/bgd-pfn/$JOB_NAME.log
#$ -j y

for i in {1..1}
do
for k in {1..1}
do
qsub umx_on_vertex.sh  "$i" "$k"
sleep 30
done
sleep 60
done
