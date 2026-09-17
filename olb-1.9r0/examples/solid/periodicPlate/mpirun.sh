#!/bin/bash

for n in 20 30 40 50 60 80 100; do
	mpirun -np 4 periodicPlate --res $n
done
