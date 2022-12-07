#!/bin/bash
#$ -N matlab_vertex_loop
#$ -o /cbica/projects/bgd-pfn/$JOB_NAME.log
#$ -j y
#n_pfn = $1
#n_vertex = $2 == 59412

for (( n=$1; n<=$2; n++))
do
  $k_start=1
  for  (( k=$3; k<=$4; k+=1000))
  do
    k_end=$k
    if [[$k_end>59412]];
    then
      k_end=59412
      qsub matlab_test.sh "$n_start" "$k_start" "$k_end"
      break
    fi
    qsub matlab_test.sh "$n_start" "$k_start" "$k_end"
    k_start=$k_end++
  done
done
