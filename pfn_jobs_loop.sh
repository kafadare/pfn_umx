#!/bin/sh
for i in {1..2}
do
for k in {1..2}
qsub umx_on_vertex.sh $i $k
sleep 30
done
sleep 60
done
