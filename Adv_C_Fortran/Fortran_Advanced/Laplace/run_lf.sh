#!/bin/bash
gfortran LaplaceEq.f90 -o lf
./lf
gnuplot -persist << EOFMarker
    set pm3d
    set hidden3d
    set size ratio 1
    splot "data" with lines
EOFMarker
