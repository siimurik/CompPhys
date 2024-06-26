#!/bin/bash

# Clean up any previous executable
rm -f benchmark_gfortran benchmark_ifort benchmark_ifx

# Compile with gfortran with optimization
echo "Compiling with gfortran..."
/usr/bin/time -v gfortran -O3 -o benchmark_gfortran benchmark.f90

# Compile with ifort with optimization
echo "Compiling with ifort..."
/usr/bin/time -v ifort -O3 -o benchmark_ifort benchmark.f90

# Compile with ifx with optimization
echo "Compiling with ifx..."
/usr/bin/time -v ifx -O3 -o benchmark_ifx benchmark.f90

# Run the gfortran executable
echo "Running gfortran executable..."
/usr/bin/time -v ./benchmark_gfortran

# Run the ifort executable
echo "Running ifort executable..."
/usr/bin/time -v ./benchmark_ifort

# Run the ifx executable
echo "Running ifx executable..."
/usr/bin/time -v ./benchmark_ifx

