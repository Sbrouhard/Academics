import numpy as np
import matplotlib.pyplot as plt
import random
import scipy

# Plots the bifurcation diagram between 2 values for a
def plot_bifurction(starting_a, final_a):
    fig, (lyaponov, bifurcation) = plt.subplots(2)

    total_iterations = 10 ** 2

    x_values = np.zeros(total_iterations + 1)
    l_values = []

    initial_condition = random.random()

    step_size = 10 ** -4

    for a in np.linspace(starting_a, final_a, 100000):
        # print(f'current a value: {a}')
        x_values[0] = initial_condition
        l_total = 0

        for i in range(0, total_iterations):
            x_values[i + 1] = a * x_values[i] * (1 - x_values[i])
            l_total += np.log(np.abs(-2 * a * x_values[i] + a)) 

        l_values.append((1.0/total_iterations) * l_total)

        longest_period = 5
        plot_orbit = longest_period
        tolerance = np.power(10.0, -5)
        inside_loop = False
        difference_between_points = np.zeros(100000000)

        for i in range(0,longest_period):
            difference_between_points[i + 1] = np.absolute(x_values[x_values.size - 1] - x_values[[x_values.size - np.power(2, longest_period)]])
            if difference_between_points[i + 1] < tolerance and not inside_loop:
                longest_period = i
                inside_loop = True

        plot_values = x_values[-1-(2**longest_period):-1]
        aVec = [a] * len(plot_values)
        bifurcation.plot(aVec, plot_values,'o', markersize=0.01,color='black')


    title = "Lyapunov Exponents iterated from " + str(x_values[0])


    x = np.linspace(starting_a, final_a, len(l_values))

    lyaponov.plot(x,l_values)
    lyaponov.set_xlabel("a parameters")
    lyaponov.set_ylabel("Lyapunov Exponent")
    lyaponov.axhline(y = 0, color = "r")
    lyaponov.title.set_text(title)

# Takes a random point in the unit interval and iterates it under ax(1-x) for a=3.9 
# For a specified number of times
def find_attracted_point_logistic(iterations) -> float:
    x = np.random.rand(1)
    for i in range(0, iterations):
        x = 3.9 * x * (1 - x)
    return x

# Given 2 points, iterates the 3.9x(1-x) 10^4 times and plots the
def find_difference_between_points(a, b):
    difference = [1]
    x = [1]
    for i in range(0, 10000):
        a = 3.9 * a * (1 - a)
        b = 3.9 * b * (1 - b)
        if(np.abs(a-b) < get_minumum(difference)):
            difference.append(float(np.abs(a - b)))
            x.append(x[len(x) - 1] + 1)
    return x[1:], difference[1:]

# Returns the minimum element of an array. 
# Numpy wasn't behaving so I just wrote a function myself
def get_minumum(x):
    min = 100000
    for element in x:
        if element < min:
            min = element
    return min


# Given a list of lists, returns the list with lowest length. Used for reconciling
# The min-distance arrays of different lengths later on
def get_array_of_min_length(x):
    min = 100000
    for array in x:
        if len(array) < min:
            min = len(array)
    return min

# Function to perform regression on. This can be ignored
def monoExp(x, m, t, b):
    return m * np.exp(-t * x) + b


# Picks a point a, compares the minimum distance plot between a and a set of 500 b's
# Which are also in the attractor. Averages them, performs regression, then plots them
def compare_points():
    pointa = find_attracted_point_logistic(1000000)

    iterations = []
    differences = []
    for i in range(0, 500):
        print(i)
        iteration, difference = find_difference_between_points(pointa, find_attracted_point_logistic(100000))
        iterations.append(iterations)
        differences.append(difference)

    # Ignoring points b which "diverge" by removing difference lists with not many elements
    new_differences = []
    for array in differences:
        if len(array) < 14:
            pass
        else:
            new_differences.append(array)
    
    differences = new_differences
            





    # Averaging the min-distances of all points
    shortest_array_length = get_array_of_min_length(differences)
    print(shortest_array_length)
    final_average = []
    for i in range(0, shortest_array_length):
        total = 0
        for array in differences:
            total += array[i]
        final_average.append(total/len(differences))

    final_iter = np.arange(0, len(final_average))
    
    # Scatter plot the points
    myplot = plt.scatter(x=final_iter, y=final_average)
    
    # Regression
    params, cv = scipy.optimize.curve_fit(monoExp, final_iter, final_average)
    m, t, b = params
    
    # Plotting
    plt.plot(final_iter, monoExp(final_iter, m, t, b), '--', label="fitted")
    # plt.yscale('log')
    plt.title(f'FUNCTION FOR RELATIONSHIP: ({m} * x^-{t}) + {b}')
    plt.show()
    


#------ UNCOMMENT FOR BIFURCATION WORK ------#
# plot_bifurction(starting_a=1, final_a=4)
# plot_bifurction(starting_a=3.90655, final_a=3.9067)
    
# ----- UNCOMMENT FOR THE MIN-DISTANCE WORK ----
compare_points()


plt.show()


