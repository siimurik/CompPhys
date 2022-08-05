echo Compilation and execution of the main Fortran file.
gfortran shoot2.f90 -o shoot2 
./shoot2
cat desc_n_calc_time.txt
echo Plotting the data with Python.
python3 graph_shoot2.py