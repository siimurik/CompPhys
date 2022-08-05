clear
    A = [   0.0    1.0    2.0
            2.0    1.0    4.0
            2.0    4.0    6.0 ]
    b = [4. 3. 7.]
    tic
    x = linsolve(A,b')
    time = toc