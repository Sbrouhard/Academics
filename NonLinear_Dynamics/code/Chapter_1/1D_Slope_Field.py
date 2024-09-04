import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

x = sp.symbols('x')   

def plot_slope_field(fun: sp.Function, x_bottom, x_top, x_resolution, end_time, time_resolution):

    #Make the function numeric
    fun_numeric = sp.lambdify(x, fun, 'numpy')

    # Get time steps
    timesteps = np.linspace(0, end_time, time_resolution)

    # Find the slope at each value
    slopes = []
    x_range = np.linspace(x_bottom, x_top, x_resolution)
    
    for value in x_range:
        slopes.append((value, fun_numeric(value)))

    plt.figure()
    
    # Build a list of lines centered at slope[0] with slope: slope[1]
    lines = []
    
    slope_length = 0
    slope_length = 0.2
    for timestep in timesteps:
        for slope in slopes:
            # Normalizng dx and dy
            angle = np.arctan(slope[1])
            dx = slope_length * np.cos(angle)
            dy = slope_length * np.sin(angle) 

            
            x0 = timestep - 0.5 * slope_length
            y0 = slope[0]
            x1 = x0 + dx
            y1 = y0 + dy

            plt.plot([x0, x1], [y0, y1], '-')
    
    plt.xlim(0, end_time)
    plt.ylim(x_bottom, x_top)
    xdotequals = r'\dot{x}='
    plt.title(f"Slope Field For ${xdotequals + sp.latex(fun)}$")




    
    


fun_to_eval = x + sp.exp(-1 * x)
plot_slope_field(fun=fun_to_eval, x_bottom=-3, x_top=3, x_resolution=30, end_time=5, time_resolution=20)
plt.show()
