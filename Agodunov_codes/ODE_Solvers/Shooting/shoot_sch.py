E_precision = 0.000001
lower_bound = 0.0
upper_bound = 4.0
E = upper_bound
dE = 1
while dE> E_precision:
    for i in lin[0:-1]:
        if i==0:
            psi[i+1]=f0[i]+dx*dpsi_0
        else: