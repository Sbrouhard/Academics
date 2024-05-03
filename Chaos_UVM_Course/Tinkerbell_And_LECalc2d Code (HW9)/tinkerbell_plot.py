import matplotlib.pyplot as plt
import numpy as np
import time
import math


def main():
    # Number of iterations of the map
    iterations = 200

    # Array containing values of the function iterated
    values = np.zeros((2, iterations), dtype=np.float64)

    # Initial conditions and placing it in values array
    initial_condition = [0.1, 0.1]
    values[0,0] = initial_condition[0]
    values[1,0] = initial_condition[1]


    plt.ion()

    figure, ax = plt.subplots(figsize=(10, 8))
    line1, = ax.plot(values[0], values[1], 'o', markersize=2)
    ax = plt.gca()
    ax.set_xlim([-0.5, 0.5])
    ax.set_ylim([-0.5, 0.5])


    c_4_perterubation = np.linspace(start= -0.5, stop=2, num=10)
    # c_4_perterubation = [0.5]

    for j in range(0, len(c_4_perterubation)):
        values = np.zeros((2, iterations))
        initial_condition = [0.1, 0.1]
        values[0,0] = initial_condition[0]
        values[1,0] = initial_condition[1]

        for i in range(0, iterations - 1):
            ax.set_ylabel("") 
            if not math.isinf(values[0, i]) and not math.isinf(values[1, i]):
                # print("All GOOD")
                next_values = tinkerbell(values[0, i], values[1, i], c_4_perterubation[j])
                print(next_values)
                values[0, i+ 1] = next_values[0]
                values[1, i+1] = next_values[1]
                # print(next_values)
                line1.set_xdata(values[0])
                line1.set_ydata(values[1])  
                ax.set_xlabel(f'C4 Value of {c_4_perterubation[j]}') 
                figure.canvas.draw()
                figure.canvas.flush_events()
            else:
                ax.set_ylabel(f'WENT TO INFINITY') 
                figure.canvas.draw()
                figure.canvas.flush_events()
                time.sleep(1)
                break

        
        
        



    

    
    


# Function Information
def tinkerbell(x_value, y_value, c_4_value):
    # Defining constant parameters of function
    c_1 = -0.3
    c_2 = -0.6
    c_3 = 2
    c_4 = c_4_value

    value_return = np.zeros((2))

    x1 = x_value**2 - y_value**2 + (c_1 * x_value) + (c_2 * y_value)
    value_return[0] = x1

    y1 = (2 * x_value * y_value) + (c_3 * x_value) + (c_4 * y_value)
    value_return[1] = y1

    return value_return

main()