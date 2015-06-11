#!/bin/bash

lieu1="$PWD"
  for i in `seq 1`
  do
    echo "iteration_${i}"
    qsub  -l ct=40:00:00 -l vmem=16G -l fsize=3G -j y -o $lieu1/out.txt $lieu1/reto.sh ${i} $lieu1
  done


