#!/bin/bash

ci="1000"

qcc -O2 -Wall findHF.c -o findHF -lm
python3 findHF.py $ci
