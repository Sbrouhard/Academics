import matplotlib.pyplot as plt
import numpy as np



def find_sign_change(domain, yvals):
    sign_changes = []

    near_0_case = False

    for i in range(1, len(domain)):
        if np.abs(yvals[i-1]) == 0 :
            sign_changes.append(domain[i])
            near_0_case = True
            
        elif (yvals[i-1] < 0 and yvals[i]) > 0 or (yvals[i-1] > 0 and yvals[i] < 0):
            sign_changes.append((domain[i-1] + domain[i])/2)
            near_0_case = True
        if np.abs(yvals[i - 1]) < 0.000000001 and not nearby_sign_change(domain[i-1], sign_changes):
            sign_changes.append(domain[i-1])
    return sign_changes

def nearby_sign_change(point, sign_changes):
    for change in sign_changes:
        if np.abs(point - change) < 0.5:
            return True
    return False

def plot_function(function, domain, additional_points):

    for element in additional_points:
        np.append(domain, element)

    outputs = []
    for element in domain:
        outputs.append(function(element))
    sign_changes = find_sign_change(domain, outputs)



    plt.scatter(domain, outputs, marker=".", s=1)
    plt.axvline(x=0, c="black")
    plt.axhline(y=0, c="black")

    arrow_points = []
    
    if sign_changes != []:
        arrow_points.append((domain[0] + sign_changes[0])/2)
        for i in range(1, len(sign_changes)):
            arrow_points.append((sign_changes[i-1] + sign_changes[i])/2)
        arrow_points.append((domain[len(domain) - 1] + sign_changes[len(sign_changes) - 1])/2)
    else:
        arrow_points.append((domain[0] + domain[int(len(domain)/2)])/2)
        arrow_points.append((domain[int(len(domain)/2)] + domain[len(domain) - 1])/2)

    arrows = []
    # print(arrow_points)
    for arrow_point in arrow_points:
        if function(arrow_point) > 0:
                plt.text(arrow_point, 0, ">", ha='center',va='center', fontsize=20)
                arrows.append(1)
        elif function(arrow_point) < 0:
                plt.text(arrow_point, 0, "<", ha='center', va='center', fontsize=20)
                arrows.append(-1)

    # print(arrows)
    print(sign_changes)
    for i in range(0, len(sign_changes)):
        try:
            if arrows[i] > 0 and arrows[i+1] < 0:
                plt.text(sign_changes[i], 0, "\u25CF", ha='center', va='center', fontsize=14)
            if arrows[i] < 0 and arrows[i+1] > 0:
                 plt.text(sign_changes[i], 0, '\u25CF', ha='center', va='center', fontsize=14, color='white', fontweight='bold')  # filled circle
                 plt.text(sign_changes[i], 0, '\u25CB', ha='center', va='center', fontsize=14, color='black')  # circle outline
            if arrows[i] > 0 and arrows[i+1] > 0:
                plt.text(sign_changes[i], 0, '\u25CF', ha='center', va='center', fontsize=14, color='white', fontweight='bold')  # filled circle
                plt.text(sign_changes[i], 0, "\u25D0", ha='center', va='center', fontsize=14)
            if arrows[i] < 0 and arrows[i+1] < 0:
                plt.text(sign_changes[i], 0, '\u25CF', ha='center', va='center', fontsize=14, color='white', fontweight='bold')  # filled circle
                plt.text(sign_changes[i], 0, "\u25D1", ha='center', va='center', fontsize=14)
        except:
             pass
    
    plt.ylim(top=y_lim_top, bottom = y_lim_bot)
    plt.title("Phase Diagram")
    

    plt.show()



x_lim_bot = -2
x_lim_top = 3
y_lim_bot = -2
y_lim_top = 2

def function(x):
    a=0
    return_value = a*x-x**3
    return return_value

plot_function(function=function, domain=(np.linspace(x_lim_bot, x_lim_top, 10000)), additional_points=[0])

