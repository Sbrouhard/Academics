import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from typing import Callable


x=sp.Symbol('x')

def euler_step(func, x0, step_size):
    slope_at_x0 = func(x0)
    return x0 + (slope_at_x0 * step_size)

def euler_estimate(func, x0, goal_time, stepsize):
    f = sp.lambdify(x, func, 'numpy')

    current_time = 0
    current_x = x0
    ts = []
    xs = []
    while current_time <= goal_time:
        current_time += stepsize
        current_x = euler_step(f, current_x, stepsize)
    return current_x

def improved_euler_step(func, x0:float, step_size:float):
    slope_x0 = func(x0)
    probe = x0 + (slope_x0 * step_size)
    estimate = x0 + ((0.5 * (slope_x0 + func(probe))) * step_size)
    return estimate

def improved_euler_estimate(func, x0, goal_time, stepsize):
    f = sp.lambdify(x, func, 'numpy')
    current_time = 0
    current_x = x0
    ts = []
    xs = []
    while current_time <= goal_time:
        current_time += stepsize
        current_x = improved_euler_step(f, current_x, stepsize)
    return current_x

def runge_kutta_step(func, x0, step_size):
    k1 = func(x0) * step_size
    k2 = func(x0 + 0.5*k1) * step_size
    k3 = func(x0 + 0.5*k2) * step_size
    k4 = func(x0 + k3) * step_size

    return x0 + (1/6) * (k1 + 2* k2 + 2 * k3 + k4)
    
def runge_kutta_fourth_order(func, x0, goal_time, stepsize):
    f = sp.lambdify(x, func, 'numpy')
    current_time = 0
    current_x = x0
    ts = []
    xs = []
    while current_time <= goal_time:
        current_time += stepsize
        current_x = runge_kutta_step(f, current_x, stepsize)
    return current_x

def approx_for_decreasing_step_size(func, x0, target_time, step_size_exp, estimation):
    t_exponents = []
    approximations = []
    print("")
    for i in range(0, step_size_exp):
        approximation = estimation(func, x0, target_time, 10**(-i))
        t_exponents.append(i)
        approximations.append(approximation)
        print(f'Step Size: 10^-{i}, Approximation: {approximation}')
    print("")
    return t_exponents, approximations
    
def plot_error_vs_t(func, x0, target_time, actual_value, step_size_exp, approx_function):
    ts, approximations = approx_for_decreasing_step_size(func, x0, target_time, step_size_exp, approx_function)
    error = []
    for approx in approximations:
        error.append(np.abs(actual_value - approx))
    plt.scatter(ts, error, marker=".", s=50)
    xdot = r'\dot{x}'
    plt.title(f"Error for Euler Method on function: ${xdot + sp.latex(func)}$")
    plt.xlabel("Scale is N where $\Delta t = 10^{-N}$")
    plt.ylabel("Error")
    plt.plot(ts, error, color='blue')

    plt.show()

def plot_error_vs_t_ln(func, x0, target_time, actual_value, step_size_exp, approx_function):
    ts, approximations = approx_for_decreasing_step_size(func, x0, target_time, step_size_exp, approx_function)
    error = []

    # Calculate log of errors
    for approx in approximations:
        error.append(np.log(np.abs(actual_value - approx)))

    # Calculate log of delta t
    newts = []
    for t in ts:
        newts.append(np.log(10**(-t)))
    
    print("")
    print(error)
    print(newts)
    print("")

    plt.scatter(newts, error, marker=".", s=50)
    xdot = r'\dot{x}'
    plt.title(f"Logarithmic Error for Euler Method on function: ${xdot + sp.latex(func)}$")
    plt.xlabel("ln(t) where $\Delta t = 10^{-n}$")
    plt.ylabel("$\ln E$")
    plt.plot(newts, error, color='blue')

    plt.show()

fun_to_test = x + sp.exp(-1 * x)
# plot_error_vs_t(func=fun_to_test, x0=1, target_time=1,actual_value=(1/np.e), step_size_exp=7, approx_function=euler_estimate)
# plot_error_vs_t_ln(func=fun_to_test, x0=1, target_time=1,actual_value=(1/np.e), step_size_exp=7, approx_function=euler_estimate)
approx_for_decreasing_step_size(fun_to_test, 0, 1, 6, runge_kutta_fourth_order)