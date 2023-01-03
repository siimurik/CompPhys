#!/bin/bash
gfortran -o invmat inverse_matrix.f90 -llapack -lblas
./invmat
gnuplot -persist << EOFMarker
    set grid
    plot 'output.csv' with lines
EOFMarker