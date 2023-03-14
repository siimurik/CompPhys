import numpy as np

val1, val2 = 0, 0
for x in range(2, 11):
    #print(x, "\t",3*x/((x**2-2)*np.log(x)))
    val1 = val1 + 3*x/((x**2-2)*np.log(x))
    val2 =  val2 + 3/(x**2*np.log(x))
    print(val1, "\t", val2)