#!/bin/bash

# Compile the Fortran program with optimization flags and linking LAPACK and BLAS
gfortran -O2 anharmonic.f90 -llapack -lblas -o an

# Iterate over values of N
for N in 4 8 12 16 24 32; do
    # Run the Fortran program with N and echo the result to the data file
    ( echo $N ; echo 0.9 ) | ./an >> data
done

# Print the selected columns from the data file using grep and awk
grep ^EV data | awk '{ print $2, $4 }'

# Start gnuplot and plot the inverse of the second column against the fourth column
gnuplot <<EOF
plot "<grep ^EV data | awk '{print 1/$2, $4}'"
EOF
