import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import time
from helper_functions import *
from matplotlib.colors import LinearSegmentedColormap



colors = [(0, 'blue'), (0.5, 'cyan'), (1, 'red')]
custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)

np.set_printoptions(formatter={'all': lambda x: "{:.4g}".format(x)})



# Given the lorenz equations, uses rk2 to take a bunch of timesteps and returns 
def iterate_lorenz_equations(initial_condition, parameters={"sigma": 10,"r": 28,"b":  8/3}, iterations=1000):

    lorenz_equations = lambda x, parameters: [-1 * parameters["sigma"] * (x[0] - x[1]), (-1 * x[0]*x[2]) + parameters["r"]* x[0] - x[1], x[0]*x[1]-parameters["b"]*x[2]]

    iterations = rk2(equations=lorenz_equations, parameters= parameters, initial_condition=initial_condition, iterations= iterations, timestep=0.01)
    
    return iterations

# Used to step the differential equations forward in time repeatedly
def rk2(equations, parameters, initial_condition, iterations, timestep):
   
    x0 = initial_condition[0]
    y0 = initial_condition[1]
    z0 = initial_condition[2]



    return_array = [[], [], []]
    for i in range(iterations):  # 10 steps as an example
        # Evaluate the slopes
        k1 = equations([x0, y0, z0], parameters)
        k1_x = k1[0]
        k1_y = k1[1]
        k1_z = k1[2]
        
        k2 = equations([x0 + 0.5 * timestep * k1_x, y0 + 0.5*timestep*k1_y, z0 + 0.5*timestep*k1_z], parameters)
        k2_x = k2[0]
        k2_y = k2[1]
        k2_z = k2[2]
        
        # Update x, y, and z using the RK2 formula
        x1 = x0 + timestep * k2_x
        y1 = y0 + timestep * k2_y
        z1 = z0 + timestep * k2_z


        # Update initial conditions for the next step
        x0 = x1
        y0 = y1
        z0 = z1
        
        # Print the values at each timestep and percent to finish
        # if i % 100 == 0:
        #     print(f"{100 * (i / iterations)} Percent Complete")
        return_array[0].append(float(x0))
        return_array[1].append(float(y0))
        return_array[2].append(float(z0))


    return return_array


# Given a set of points, computes their boxcounting dimension
def compute_boxcounting_dimension(points):
    
    minimum_x = np.min(points[0])
    minimum_y = np.min(points[1])
    minimum_z = np.min(points[2])

    max_x = np.max(points[0])
    max_y = np.max(points[1])
    max_z = np.max(points[2])

    new_array = []
    for i in range(0, len(points[0])):
        new_array.append([points[0][i], points[1][i], points[2][i]])
    points = new_array

    total_width = max_x - minimum_x
    total_length = max_y - minimum_y
    total_height = max_z - minimum_z

   



    maximum_step = 5
    x = np.zeros(maximum_step + 1)
    y = np.zeros(maximum_step + 1)

    for step in range(1, maximum_step + 1):
        print(f'Current subdivision number: {step}')
        number_of_boxes = 0
        number_of_subdivisions = 2 ** step

        box_width = total_width / number_of_subdivisions
        box_height = total_height / number_of_subdivisions
        box_length = total_length / number_of_subdivisions

        for x_subdivision in range(1, number_of_subdivisions + 1):
            lower_left_box_x = minimum_x + ((x_subdivision - 1) * box_width)
            upper_right_box_x = minimum_x + ((x_subdivision) * box_width)
            indicies_of_xs_in_box = find_indicies(points, lambda x: x[0] >= lower_left_box_x and x[0] < upper_right_box_x)
            points_filtered_x = list(filter(lambda x: points.index(x) in indicies_of_xs_in_box, points))
            
            for y_subdivision in range(1, number_of_subdivisions + 1):
                lower_left_box_y = minimum_y + ((y_subdivision - 1) * box_length)
                upper_right_box_y = minimum_y + ((y_subdivision) * box_length)
                
                indicies_of_ys_in_box = find_indicies(points_filtered_x,lambda x: (x[1] >= lower_left_box_y) and (x[1] < upper_right_box_y))
                points_filtered_xy = list(filter(lambda x: points_filtered_x.index(x) in indicies_of_ys_in_box, points_filtered_x))
                
                if len(points_filtered_xy) > 0:
                    number_of_boxes += 1
                for z_subdivision in range(1, number_of_subdivisions + 1):
                    lower_left_box_z = minimum_z + ((z_subdivision - 1) * box_height)
                    upper_right_box_z = minimum_z + (z_subdivision * box_height)
                    indicies_of_zs_in_box = list(find_indicies(points_filtered_xy,lambda x: (x[2] >= lower_left_box_z) and (x[2] < upper_right_box_z)))
                    points_filtered_xyz = list(filter(lambda x: points_filtered_xy.index(x) in indicies_of_zs_in_box, points_filtered_xy))
                    if len(points_filtered_xyz) > 0:
                        number_of_boxes += 1
        x[step] = step * np.log(2)
        y[step] = np.log(number_of_boxes)

  



    m,b = np.polyfit(x, y, 1)
    # Returning the slope of the line of best fit for the bcd. This is the dimension!

    plt.scatter(x[1:], y[1:])
    plt.title("Box Counting Dimension Estimate")

    return m

    

# Computes the lyapunov exponent of the lorenz attractor given an initial condition and
# Dictionary containing the parameters
def lyapunov_exponent_lorenz(v0, parameters, number_iterations=1000):    
    # Dimension of the initial condition, iterations, base state of exponent, timestep, and setting v to be our IC
    m = len(v0)
    number_iterations = number_iterations
    lyapunov_exponent = np.zeros(m)
    timeStep = 0.01
    v = v0

    # Defining the lorenz equations and its jacobian as anonymous functions. they can be called given an IC x and a dictionary of parameters
    lorenz_equations = lambda x, parameters: np.array([-1 * parameters["sigma"] * (x[0] - x[1]), (-1 * x[0]*x[2]) + parameters["r"]* x[0] - x[1], x[0]*x[1]-parameters["b"]*x[2]])
    lorenz_jacobian = lambda x, parameters: np.matrix([[-1 * parameters["sigma"], parameters["sigma"], 0], [(-1 * x[2]) + parameters["r"], -1, -1 * x[0]], [x[1], x[0], -1 * parameters["b"]]])
    
    # This is just a function which takes 2 arrays and multiplies them. Will be used to deform the unit ball with the jacobian
    Jt = lambda x, A: A @ x
    q = np.eye(m)

    # The fun stuff. Morph the unit ball until we get the deformed unit ball 1 time step in the future
    # Then do QR decomp on it, get the eigenvalues of the upper-triangular, and calc the
    # Lyapunov from there :)
    for i in range(1, number_iterations):
        DF1 = np.eye(m)
        for j in range(1, int(1/timeStep)):
            DF1 = stepIt(Jt, DF1, lorenz_jacobian(v, parameters), timeStep)
            v = stepIt(lorenz_equations, v, parameters, timeStep)

        Z = DF1 @ q
        [q, r] = np.linalg.qr(Z, mode='complete')
        lyapunov_exponent += np.log(np.abs(r.diagonal()))/number_iterations
    
    return lyapunov_exponent
    

# This takes a function, condition, parameters, and a timestep and "shifts" the IC forward by the amount
# of time specified in timestep. 
def stepIt(function, x, p, timeStep):
    s1 = np.array(function(x, p))
    s2 = np.array(function(x + timeStep * (s1/2), p))
    s3 = np.array(function(x + timeStep * (s2/2), p))
    s4 = np.array(function(x + timeStep * s3, p))
    x = np.array(x + (timeStep*(s1 + 2*s2 + 2*s3 + s4)/6))
    return x


# Computes the max lyapunov exponent along all values of sigma and b. mode is either "read"
# or "write" depending if you already have a file with data
def plot_max_exponent_sigma_b_plane(mode):
    
    # If we want to generate new points: 
    if mode == "write":
        write_file = open("sigma_b_plane_exponents.csv", "w")
        exponent_calculation = []
        sigmas = np.linspace(0, 60, 120)
        bs = np.linspace(0, 10, 50)

        # Iterate through sigma and b, calculating the exponent and writing it to file
        for sigma in sigmas:
            print(f'{sigma/60 * 100} complete')
            for b in bs:
                params = {"sigma": sigma, "r": 28, "b": b}
                v0 = [1.508870, -1.531271, 25.46091]

                max_exponent = np.max(lyapunov_exponent_lorenz(v0, params, number_iterations=100))
                exponent_calculation.append([sigma, b, max_exponent])
                write_file.write(f"{sigma}, {b}, {max_exponent}\n")
        write_file.close()
    
    # Takes the data from file and converts to a numpy array
    data = np.genfromtxt('sigma_b_plane_exponents.csv', delimiter=',', dtype=[('x', float), ('y', float), ('z', float)])
    array_data = np.array([(row['x'], row['y'], row['z']) for row in data])

    # Plotting the data with red points being chaotic
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(array_data[:, 0], array_data[:, 1], array_data[:, 2], marker=".", s=10, c=color_over_0(array_data[:, 2]))
    plt.show()


# Given an orbit, finds the local max values and plots them in a time-series fashion
def plot_successive_z_maxima(orbit):
    orbit_array = np.array(orbit)
    z_values = orbit_array[2, :]

    z_max = get_local_maxima(z_values)
    z_max_offset = z_max[1:]


    plt.scatter(z_max[:len(z_max)- 1], z_max_offset, marker=".", s= 2)
    plt.show()
    
# Iterates r and keeps track of non-transient orbit. Plots the points which remain!
def plot_bifurcation_of_z_max():
    r_values = np.linspace(0, 325, 10000)
    points_to_plot = []

    for r in r_values:
        parameters = {"sigma": 10, "r": r, "b": (8/3)}
        v0 = [1.508870, -1.531271, 25.46091]
        orbit = iterate_lorenz_equations(v0, parameters, iterations=1000)

        orbit_zs = np.array(orbit)[2, :]
        z_maxs = get_local_maxima(orbit_zs[500:])

        for z_max in z_maxs:
            points_to_plot.append([r, z_max])
    
    points_to_plot = np.array(points_to_plot)
    plt.scatter(points_to_plot[:, 0], points_to_plot[:, 1], marker=".", s= 0.001)
    plt.xlim(0, 325)
    plt.ylim(0, 400)
    plt.show()


def generate_heat_map(r_value, z_value):
    lorenz_jacobian = lambda x, parameters: np.matrix([[-1 * parameters["sigma"], parameters["sigma"], 0], [(-1 * x[2]) + parameters["r"], -1, -1 * x[0]], [x[1], x[0], -1 * parameters["b"]]])
    
    parameters = {"sigma": 10, "r": r_value, "b": (8/3)}
    xs = np.linspace(-10, 10, 100)
    ys = np.linspace(-10, 10, 100)

    points_to_plot = []

    for x in xs:
        for y in ys:
            jacobian = lorenz_jacobian([x, y, z_value], parameters)
            eigenvalues = np.linalg.eigvals(jacobian)
            real_components = np.real(eigenvalues)
            max_real_eig = np.max(real_components)

            points_to_plot.append([x, y, max_real_eig])
    points_to_plot = np.array(points_to_plot)
    
    xs = points_to_plot[:, 0]
    ys = points_to_plot[:, 1]
    zs = points_to_plot[:, 2]
    return[xs, ys, zs]

def plot_heat_map(r_val):
    fig, axs = plt.subplots(3, figsize=(4, 8), gridspec_kw={'hspace': 0.4})
    z_15 = generate_heat_map(z_value=15, r_value=r_val)
    z_25 = generate_heat_map(z_value=25, r_value=r_val)
    z_35 = generate_heat_map(z_value=35, r_value=r_val)

    fig1 = axs[0].scatter(z_15[0], z_15[1], c = z_15[2], cmap="viridis", vmin = -3, vmax = 5)
    axs[0].set_title(f'Z = 15, r = {r_val}', fontsize=8)


    fig2 = axs[1].scatter(z_25[0], z_25[1], c = z_25[2], cmap="viridis",  vmin = -3, vmax = 5)
    axs[1].set_title(f'Z = 25, r = {r_val}', fontsize=8)

    fig3 = axs[2].scatter(z_35[0], z_35[1], c = z_35[2], cmap="viridis",  vmin = -3, vmax = 5)
    axs[2].set_title(f'Z = 35, r = {r_val}', fontsize=8)


    # Add colorbars
    cbar = fig.colorbar(fig3, ax=axs, orientation='horizontal')
    cbar.set_label('Real Component of EVs', fontsize=8)
        
# ------- Homework Solutions! Uncomment the code you're interested in -----------#
def main():

    # -----------  Question 1: Compute the fractal dimension of the lorenz attractor -------------------------------#
    print("ESTIMATING DIMENSION")
    v0 = [1.508870, -1.531271, 25.46091]
    parameters={"sigma": 10, "r": 24, "b": (8/3)}

    orbit = iterate_lorenz_equations(v0, parameters, iterations=10000)
    box_counting_dimension = compute_boxcounting_dimension(orbit)
    print(f'\n\nBOX COUNTING DIMENSION: {box_counting_dimension}\n\n')
    plt.xlabel(f"Estimated Dimension: {box_counting_dimension}")
    plt.show()

    input("Press Enter to continue to Question 2")


    # ----------- Question 2: Calculate the Lyapunov Exponents for the Lorenz Equation with r=12, 24.5, 28 ----------#
    print("Calculating Lyapunov Exponents for varying R values. Sigma: 10, B: 8/3")

    v0 = [1.508870, -1.531271, 25.46091]
    exp_r_12 = lyapunov_exponent_lorenz(v0, {"sigma": 10, "r": 12, "b": (8/3)})
    print(f'R = 12:\n Exponents: {exp_r_12} \n')

    exp_r_24_5 = lyapunov_exponent_lorenz(v0, {"sigma": 10, "r": 24.5, "b": (8/3)})
    print(f'R = 24.5:\n Exponents: {exp_r_24_5} \n')

    exp_r_28 = lyapunov_exponent_lorenz(v0, {"sigma": 10, "r": 28, "b": (8/3)})
    print(f'R = 28:\n Exponents: {exp_r_28} \n')

    input("\nPress Enter to continue to Question 2\n\n")

    # ----------- Question 3: Plot the max lyapunov exponent as you pass through a plane of values 
    print("Plotting maximum lyapunov exponent on plane of sigma and b. Close the figure to move to Question 4")
    plot_max_exponent_sigma_b_plane(mode="read")

    input("Press Enter to Continue to Question 4")

    # ----------- Question 4: Plot the successive z_maxima and bifurcation diagram
    v0 = [1.508870, -1.531271, 25.46091]
    parameters={"sigma": 10, "r": 25, "b": (8/3)}

    orbit = np.array(iterate_lorenz_equations(v0, parameters, iterations=100000))
    plot_successive_z_maxima(orbit)

    print("Plotting bifurcation diagram now...")
    plot_bifurcation_of_z_max()

    input("Press Enter to Continue to Question 5")

    # ---------- Quesiton 5: Meteorology! 

    plot_heat_map(r_val=12)
    plot_heat_map(r_val=24)
    plot_heat_map(r_val=35)  
    plt.show()



main()