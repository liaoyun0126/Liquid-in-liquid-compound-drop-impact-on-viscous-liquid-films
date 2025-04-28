#!/bin/bash

MAXlevel="12"
tmax="8.0"
We="4.0"
Ohd="0.034"
Bo="0.5"
Ohf="0.670"
hf="0.35"

mkdir VideoBas2
#cd VideoBas
#rm *.png
#cd ..
qcc -O2 -Wall getData_InsideDrop.c -o getData_InsideDrop -lm
qcc -O2 -Wall getFacet1.c -o getFacet1 -lm
qcc -O2 -Wall getFacet2.c -o getFacet2 -lm
qcc -O2 -Wall getFacet3.c -o getFacet3 -lm
python3 DropImpactFilms2.py $We $Ohd $Ohf $hf
