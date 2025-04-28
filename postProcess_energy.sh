#!/bin/bash

ci="1000"
Rhor="0.001"
Ohd="0.75"
Ohf="0.670"
Bo="0.5"


qcc -O2 -Wall getEnergyAxi.c -o getEnergyAxi -lm
python3 getEnergyScript.py $ci $Rhor $Ohd $Ohf $Bo



