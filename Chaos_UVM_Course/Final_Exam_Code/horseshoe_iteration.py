import numpy as np
import matplotlib.pyplot as plt

def upper_y_horseshoe(coordinates):
    return (1 - coordinates[0]/3, 3- (3 * coordinates[1]))

def lower_y_horseshoe(coordinates):
    return (coordinates[0]/3, 3 * coordinates[1])

def upper_y_horseshoe_inverse(coordinates):
    return (-3 * coordinates[0] + 3, -1/3 * coordinates[1] + 1)

def lower_y_horseshoe_inverse(coordinates):
    return (3 * coordinates[0], 1/3 * coordinates[1])
# def lower_y_horseshoe(coordinates): 
# fig = plt.figure()
# ax = fig.subplots(2)
# plt.scatter(initials[0], initials[1])
# plt.scatter(next[0], next[1])
# plt.show()
