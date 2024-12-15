#!/bin/bash

for i in $(seq 0 10) ; do
  echo $i
	S=$((2**$i))
	make SAMPLES=32 STEP_SIZE=0.25 STEPS=$S LENGTH=300
	mv simul.h5 simul.h5-${i}
done

cp simul.h5.orig simul.h5 
