#!/bin/bash
DELTA=100000
for i in `seq 1 50`;
do
		c=$(($i * $DELTA))
                ./rgen "${c}" ~/meshes/random/"${c}".off

done  
