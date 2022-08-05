echo Compilation and execution of the main Fortran file.
gfortran Ode12.f90 -o Ode12.exe
./Ode12.exe
cat loops_n_time12.txt
echo Plotting the data with Python.
python3 graph_Ode12.py 
