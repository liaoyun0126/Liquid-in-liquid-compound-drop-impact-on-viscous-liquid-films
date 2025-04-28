#!/bin/bash

MAXlevel="12"
tmax="4.0"
We="4.0"
Ohd="0.75"
Bo="0.5"
Ohf="0.670"
hf="0.001"

#export OMP_NUM_THREADS=16
#qcc -fopenmp -Wall -O2 dropFilm.c -o dropFilm -lm
#./dropFilm $MAXlevel $tmax $We $Ohd $Bo $Ohf $hf

CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 dropFilm.c -o dropFilm -lm
mpirun -np 32 dropFilm $MAXlevel $tmax $We $Ohd $Bo $Ohf $hf
