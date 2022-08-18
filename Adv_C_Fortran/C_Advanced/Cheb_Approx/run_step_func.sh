echo Compiling the main code.
gcc -Wall -I/home/siim/gsl/include -c step_func_Cheb.c 
gcc -L/home/siim/gsl/lib step_func_Cheb.o -lgsl -lgslcblas -lm
export LD_LIBRARY_PATH="/home/siim/gsl/lib"
echo Executing and printing the results into the 'data.dat' file.
./a.out > data.dat
echo Plotting the results in Python.
python3 plot_step_func.py
