import numpy as np
import matplotlib.pyplot as plt
import random

fig, (lyaponov, bifurcation) = plt.subplots(2)


total_iterations = 10 ** 2

x_values = np.zeros(total_iterations + 1)
l_values = []

initial_condition = random.random()

starting_a_value = 1
final_a_value = 4
step_size = 10 ** -4

for a in np.arange(starting_a_value, final_a_value, step_size):
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


x = np.linspace(starting_a_value, final_a_value, len(l_values))

lyaponov.plot(x,l_values)
lyaponov.set_xlabel("a parameters")
lyaponov.set_ylabel("Lyapunov Exponent")
lyaponov.axhline(y = 0, color = "r")
lyaponov.title.set_text(title)
plt.show()