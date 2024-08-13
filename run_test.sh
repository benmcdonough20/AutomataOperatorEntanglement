#!/bin/bash

for n in {6..10..2};
do
	for h in $(seq 1 $n);
	do
		echo $n $h
	done
done
