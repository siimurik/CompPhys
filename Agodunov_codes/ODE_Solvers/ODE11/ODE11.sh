echo Compilation and execution of the main Fortran file.
gfortran Ode11.f90 -o Ode11.exe
./Ode11.exe
cat loops_n_time1.txt
echo Plotting the data with Python.
python3 graph_Ode1.py 