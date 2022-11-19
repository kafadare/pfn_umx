#!/bin/bash
#$ -N pfn_loop
#$ -o /cbica/projects/bgd-pfn/$JOB_NAME.log
#$ -j y
#n_pfn = $1
#n_vertex = $2 == 59412
for n in {1..1}
do
for k in {1..20}
do
qsub matlab_test.sh "$n" "$k"
file="/cbica/projects/bgd-pfn/vertex_columns/PFN"
file+="$n"
file+="V$k.csv"
while [ ! -f $file ]; # true if /your/file does not exist
do
  sleep 2;
done
sleep 1;
qsub umx_on_vertex.sh "$n" "$k"
sleep 30
done
sleep 60
done
