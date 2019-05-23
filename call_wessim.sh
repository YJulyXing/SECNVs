#!/bin/bash

genome=$1
region=$2
nreads=$3
read_length=$4
frag_size=$5
stdev=$6
model=$7
output=$8
qual=$9
pp=${10}

if [ $pp == "p" ]
then
	python Wessim1.py -R $genome -B $region -p -n $nreads \
	-l $read_length -f $frag_size -M $model -d $stdev \
	-o $output -t 10 -q $qual
else
	python Wessim1.py -R $genome -B $region -n $nreads \
	-l $read_length -f $frag_size -M $model -d $stdev \
	-o $output -t 10 -q $qual
fi