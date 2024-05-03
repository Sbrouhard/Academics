import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import time

# This colors the points over 0 red and the points below 0 blue
def color_over_0(z_coordinates):
    colors = []
    for z_point in z_coordinates:
        if z_point < 0:
            colors.append("blue")
        else:
            colors.append("red")
    return colors



# Plots a set of points in 3d
def plot_points(orbit, markersize = 0.2):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(orbit[0], orbit[1], orbit[2], marker=".", s=markersize)
    plt.show()


# Given a list and a function which returns a boolean, 
# returns the indicies of elements which meet the condition
def find_indicies(lst, condition):
    return_array = []
    for i in range(0, len(lst)):
        if condition(lst[i]):
            return_array.append(i)
    return return_array


# Given a list, find the local maxima
def get_local_maxima(values):
    maximal_values = []
    for i in range(1, len(values) - 1):
        if values[i - 1] < values[i] and values[i + 1] < values[i]:
            maximal_values.append(values[i])

    return maximal_values