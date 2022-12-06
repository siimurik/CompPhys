echo Compiling the main code.
gcc -I/usr/include/gsl -lgsl -lgslcblas -lm step_func_Cheb.c -o step
echo Executing and printing the results into the 'data.dat' file.
./step > data.dat
echo Plotting the results in Python.
python3 plot_step_func.py
