#!/bin/bash
REPEATS=1
DIV=$(($REPEATS*($REPEATS-1)))
prefix="noise"

DELTA=100000
for i in `seq 1 50`;
do
		c=$(($i * $DELTA))
		FILE=${prefix}""${c}".off"
		printf '%s      with    n=%i\n' ${FILE} $c
		ELAPSED=$(eval "./mdt ~/meshes/noise/${FILE}")
done
