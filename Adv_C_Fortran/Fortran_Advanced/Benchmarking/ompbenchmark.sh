#!/bin/bash

# Number of threads for OpenMP
export OMP_NUM_THREADS=4

# Paths to the compilers (adjust if necessary)
GFORTRAN=gfortran
IFORT=ifort
IFX=ifx

# Optimization flags
GFORTRAN_FLAGS="-O3 -fopenmp -o ompbenchmark_gfortran"
IFORT_FLAGS="-O3 -xHost -ipo -fast -qopenmp -diag-disable=10448 -o ompbenchmark_ifort"
IFX_FLAGS="-O3 -xHost -ipo -fast -qopenmp -o ompbenchmark_ifx"

# Compiling with gfortran
echo "Compiling with gfortran..."
/usr/bin/time -v $GFORTRAN $GFORTRAN_FLAGS ompbenchmark.f90

# Compiling with ifort
echo "Compiling with ifort..."
/usr/bin/time -v $IFORT $IFORT_FLAGS ompbenchmark.f90

# Compiling with ifx
echo "Compiling with ifx..."
/usr/bin/time -v $IFX $IFX_FLAGS ompbenchmark.f90

# Running gfortran executable
echo "Running gfortran executable..."
/usr/bin/time -v ./ompbenchmark_gfortran

# Running ifort executable
echo "Running ifort executable..."
/usr/bin/time -v ./ompbenchmark_ifort

# Running ifx executable
echo "Running ifx executable..."
/usr/bin/time -v ./ompbenchmark_ifx

