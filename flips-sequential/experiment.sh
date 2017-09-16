#!/bin/bash
DELTA=100000
for i in `seq 1 50`;
do
		c=$(($i * $DELTA))
		echo "**** performance of lawson on cgal${c}.off"
                ./app ~/random_meshes/cgal/cgal"${c}".off
		echo "***************** end ********************"
done
