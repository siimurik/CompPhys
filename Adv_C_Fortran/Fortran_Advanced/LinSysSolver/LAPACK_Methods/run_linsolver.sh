#!/bin/bash
gfortran -o dgesv dgesv_solve.f90  -llapack -lblas
./dgesv
gnuplot -persist << EOFMarker
    set grid
    plot 'out.csv' with lines
EOFMarker