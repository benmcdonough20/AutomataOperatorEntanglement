#!/bin/bash

for h in $(seq 1 5);
do
	julia finite_size_scratch.jl 10 $h
done
