echo Compiling the main code.
gcc mixing.c -o mix -lm -lfftw3
echo Executing
./mix
echo Plotting with Python
python3 plot_mixF.py 