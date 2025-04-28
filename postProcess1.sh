#!/bin/bash

MAXlevel="12"
tmax="4.0"
We="4.0"
Ohd="0.75"
Bo="0.5"
Ohf="0.670"
hf="0.001"

qcc -O2 -Wall getData_InsideDrop.c -o getData_InsideDrop -lm
qcc -O2 -Wall getFacet1.c -o getFacet1 -lm
qcc -O2 -Wall getFacet2.c -o getFacet2 -lm
qcc -O2 -Wall getFacet3.c -o getFacet3 -lm
python3 DropImpactFilms1.py $We $Ohd $Ohf $hf
