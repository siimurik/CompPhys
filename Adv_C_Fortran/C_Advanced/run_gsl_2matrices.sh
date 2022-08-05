 gcc -Wall -I/home/siim/gsl/include -c gsl_2matrices.c
 gcc -L/home/siim/gsl/lib gsl_2matrices.o -lgsl -lgslcblas -lm
 export LD_LIBRARY_PATH="/home/siim/gsl/lib"
 ./a.out
