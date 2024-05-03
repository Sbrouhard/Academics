import numpy as np
import sympy as sp
from sympy import pprint
import matplotlib.pyplot as plt



# Defining the functions we're interested in 



x = sp.Symbol('x')
y = sp.Symbol('y')
a = sp.Symbol('a')
b = sp.Symbol('b')
xy = sp.Matrix([x, y])

# Setting up the Henon Map
a = 1.4
b = 0.3
henon = sp.Matrix([a - x**2 + (b * y),  x])
henon_jac = henon.jacobian(xy)


#Setting up the tinkerbell map with arbitrary constants. Substitute in later
c1 = sp.Symbol('c1')
c2 = sp.Symbol('c2')
c3 = sp.Symbol('c3')
c4 = sp.Symbol('c4')

tinkerbell = sp.Matrix([(x**2 - y**2) + (c1 * x) + (c2 * y), (2 * x * y) + (c3 * x) + (c4 * y)])
tinkerbell_jac = tinkerbell.jacobian(xy)

# Setting up ikeda map
R = 1
C1 = 0.4
C2 = 0.9
C3 = 6
t = C1 - C3/(1 + x**2 + y**2)
ikeda = sp.Matrix([R + C2 * (x * sp.cos(t) - y * sp.sin(t)), C2 * (x * sp.sin(t) + y * sp.cos(t))])
ikeda_jac = ikeda.jacobian(xy)




def compute_exponents(function_name, function, initial_conditions):
    initial_conditions = initial_conditions
    iterations = 500
    values = np.zeros((2, iterations ), dtype = np.float64)
    values[0,0] = initial_conditions[0]
    values[1,0] = initial_conditions[1]


    Ls = np.zeros((2, iterations))
    identity_matrix = sp.eye(2)
    for i in range (0, iterations - 1):
        next_values = function.subs({x: values[0,i], y: values[1, i]})
        values[0, i+1] = next_values[0]
        values[1, i+1] = next_values[1]

        to_decomp = function.jacobian(xy).subs({x: values[0, i], y:values[1, i]}) * identity_matrix
        [identity_matrix, R] = np.linalg.qr(to_decomp, mode='complete')
        # print("R MATRIX:")
        # pprint(R)
        Ls[0, i] = R[0, 0]
        Ls[1, i] = R[1, 1]




    xs = np.log(np.abs(Ls[0][:iterations-1]))
    ys = np.log(np.abs(Ls[1][:iterations-1]))
    print(f"Lyapunov Exponents for {function_name}")
    print(np.average(xs))
    print(np.average(ys))





print()
compute_exponents("Tinkerbell with Quasi-Periodic Orbit", tinkerbell.subs({c1: -0.3, c2: -0.6, c3: 2, c4: 0.5}), initial_conditions=[0.1,0.1])
print()
compute_exponents("Henon Map", henon, initial_conditions=[0,0])
print()
compute_exponents("Tinkerbell with Chaotic", tinkerbell.subs({c1: 0.9, c2: -0.6, c3: 2, c4: 0.5}), initial_conditions=[0.1,0.1])
print()
compute_exponents("Ikeda Map", ikeda, initial_conditions=[0.1, 0.1])
print()