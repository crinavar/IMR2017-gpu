#!/bin/bash
REPEATS=50
DIV=$(($REPEATS*($REPEATS-1)))
prefix="noise"

DELTA=100000
for i in `seq 1 50`;
do
		c=$(($i * $DELTA))
		FILE=${prefix}""${c}".off"
		
		ACCUM_TIME=$((0))
		ACCUM_SQRT=$((0))
		printf '%s      with    n=%i' ${FILE} $c
		for i in `seq 1 $REPEATS`
		do
			ELAPSED=$(eval "./mdt ~/meshes/noise/${FILE}")
			#printf 'run %i time %f\n' $i $ELAPSED
			MODULUS=$(($i % 10))
			if test $MODULUS -eq 0
			then
				printf '.'
			fi
			ACCUM_TIME=`echo "scale=16; ($ACCUM_TIME + $ELAPSED)" | bc`
			ACCUM_SQRT=`echo "scale=16; ($ACCUM_SQRT + ($ELAPSED*$ELAPSED)/($REPEATS -1))" | bc`
		done
		MEAN=`echo "scale=32; ($ACCUM_TIME / $REPEATS)" | bc`
		STDEV=`echo "scale=32; sqrt($ACCUM_SQRT - ($ACCUM_TIME * $ACCUM_TIME)/$DIV)" | bc -l`
		STDERR=`echo "scale=32; ($STDEV / $MEAN)" | bc`
		printf 'mean=%f  stderr=%f\n' $MEAN $STDERR
		printf '%s   %i    %f    %f    %f\n' $FILE $c $MEAN $STDEV $STDERR >> noiseperf_mdt.dat
done
