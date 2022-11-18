#!/bin/bash
#$ -N pfn_loop
#$ -o /cbica/projects/bgd-pfn/$JOB_NAME.log
#$ -j y

for i in {1..1}
do
for k in {1..1}
do
matlab -r "get_vertex_col($i,$k);quit"
qsub umx_on_vertex.sh  "$i" "$k"
file="V"
file+="$i"
file+="_PFN$2.csv"
rm $file
sleep 30
done
sleep 60
done
