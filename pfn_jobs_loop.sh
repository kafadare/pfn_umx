#!/bin/bash
for i in {1..2}
do
for k in {1..2}
qsub umx_on_vertex.sh -v "network=$i,vertex=$k"
sleep 30
done
sleep 60
done
