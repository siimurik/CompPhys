#!/bin/bash
#echo Compiling the main code.
#gcc mixRevamp.c -o mixR -lm -lfftw3
echo Executing
./mixR
echo Plotting with Python
python3 plot_mixC.py 