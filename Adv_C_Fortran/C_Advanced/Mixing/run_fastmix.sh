#!/bin/bash

echo Compiling the main code.
gcc mixFast.c -o mixF -lm
echo Executing
./mixF
echo Plotting with Python
python3 plot_mixC.py 
rm mixF coloursF.csv