# Solve the quadratic equation ax**2 + bx + c = 0

# import complex math module
import cmath

# Taking 3 inputs
a, b, c = [float(x) for x in input("Enter three values\n").split(', ')]
print("\nThe value of a is {}, b is {} and c is {}".format(a, b, c))

# calculate the discriminant
d = (b**2) - (4*a*c)

# find two solutions
sol1 = (-b-cmath.sqrt(d))/(2*a)
sol2 = (-b+cmath.sqrt(d))/(2*a)

print('\nThe solutions are {0} and {1}'.format(sol1,sol2))
