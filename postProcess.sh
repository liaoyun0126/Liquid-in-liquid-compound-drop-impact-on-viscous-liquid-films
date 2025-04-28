#!/bin/bash

MAXlevel="11"
tmax="2.5"
We="4.0"
Ohd="0.034"
Bo="0.5"
Ohf="0.670"
hf="0.35"

mkdir VideoBas
#cd VideoBas
#rm *.png
#cd ..
qcc -O2 -Wall getData.c -o getData -lm
qcc -O2 -Wall getFacet1.c -o getFacet1 -lm
qcc -O2 -Wall getFacet2.c -o getFacet2 -lm
python3 video_Basilisk.py $We $Ohd $Ohf $hf
