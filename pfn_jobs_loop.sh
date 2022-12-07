#!/bin/bash
#$ -N pfn_loop
#$ -o /cbica/projects/bgd-pfn/$JOB_NAME.log
#$ -j y
#n_pfn = $1
#n_vertex = $2 == 59412

# for (( n=$1; n<=$2; n++))
# do
#   k_start = 1
#   for  (( k=$3; k<=$4; k+=900))
#   do
#     k_end = $k
#     qsub matlab_test.sh "$n_start" "$k_start" "$k_end"
#     k_start = $k_end + 1
#   done
# done
# file="/cbica/projects/bgd-pfn/vertex_columns/PFN"
# file+="$n"
# file+="V$k.csv"
# while [ ! -f $file ]; # true if /your/file does not exist
# do
#   sleep 2;
# done
for ((n=$1; n<=$2; n++))
do
  k=$3
  qsub umx_on_vertex.sh "$n" "$k"
  sleep 2h;
done
