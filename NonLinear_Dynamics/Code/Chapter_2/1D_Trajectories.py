import numpy as np
import matplotlib.pyplot as plt



def plot_trajectories(function, start_domain, end_domain, num_values, time_passed):
    outputs = []
    t_values = []
    initial_x_values = np.linspace(start_domain, end_domain, num_values)
    for x_val in initial_x_values:
        ts, y_values = rk2_solver(function, x_val, 0, time_passed, 0.0001)
        for i in range(0, len(ts)):
            t_values.append(ts[i])
            outputs.append(y_values[i])
    plt.scatter(t_values, outputs, 0.1, marker=".")
    plt.ylim(bottom=y_lim_bot, top = y_lim_top)
    plt.title("Trajectories of x_0's")
    plt.show()
        



def rk2_step(func, x, t, dt):
    """
    Perform a single step of the RK2 method for the ODE dx/dt = func(x, t).
    
    Parameters:
    - func: Function representing the ODE, dx/dt = func(x, t)
    - x: Current value of the solution
    - t: Current time
    - dt: Time step size
    
    Returns:
    - x_next: Estimated value of the solution at time t + dt
    """
    k1 = dt * func(x)
    k2 = dt * func(x + 0.5 * k1)
    x_next = x + k2
    return x_next

def rk2_solver(func, x0, t0, t_end, dt):
    """
    Solve the ODE dx/dt = func(x, t) using the RK2 method.
    
    Parameters:
    - func: Function representing the ODE, dx/dt = func(x, t)
    - x0: Initial value of the solution
    - t0: Initial time
    - t_end: End time
    - dt: Time step size
    
    Returns:
    - ts: Arrax of time points
    - xs: Arrax of solution values at corresponding time points
    """
    ts = [t0]
    xs = [x0]
    
    t = t0
    x = x0
    
    while t < t_end:
        x = rk2_step(func, x, t, dt)
        t += dt
        
        ts.append(t)
        xs.append(x)
    
    return ts, xs


def fun1_to_test(x):
    c=0.9
    return_value = -x
    return return_value


x_lim_bot = 1
x_lim_top = 1
y_lim_bot = -1
y_lim_top = 2
plot_trajectories(function=fun1_to_test, start_domain=x_lim_bot, end_domain=x_lim_top, num_values=1, time_passed=5)
