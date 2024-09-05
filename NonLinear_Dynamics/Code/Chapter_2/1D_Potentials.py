import sympy as sp
import numpy as np
import matplotlib.pylab as plt

x = sp.symbols('x')   


def plot_potentials(function_to_eval: sp.Function, x_bottom, x_top, y_bottom, y_top):
    # Find potential 
    potential_function = -1 * sp.integrate(function_to_eval, x)
    print(potential_function)
    interval = sp.Interval(x_bottom, x_top)

    # Find the mins and maxes
    f_prime = sp.diff(potential_function, x)
    f_double_prime = sp.diff(f_prime, x)

    zeroes = sp.solveset(f_prime, x, interval)
    maxes = []
    mins = []
    for zero in zeroes:
        if zero.is_real:
            if f_double_prime.subs(x, zero) < 0:
                maxes.append(zero)
            else:
                mins.append(zero)
    
    #turn sympy function numeric
    function_numeric = sp.lambdify(x, potential_function, 'numpy')

    #Get all points
    domain = np.linspace(x_bottom, x_top, 10000)
    range = []
    for number in domain:
        range.append(function_numeric(number))

    # Annotate mins and maxes on plot
    for max in maxes:
        max = max.evalf()
        if max.is_real:
            xpoint = float(max)
            ypoint = function_numeric(xpoint)
            plt.text(xpoint, ypoint, '\u25CF', ha='center', va='center', fontsize=14, color='white', fontweight='bold')  # filled circle
            plt.text(xpoint, ypoint, '\u25CB', ha='center', va='center', fontsize=14, color='black')  # circle outline
            # plt.annotate(f'{xpoint:0.2f}', xy=(xpoint, ypoint), xytext=(xpoint - 0.15, ypoint + 0.2), fontsize=13)


    for min in mins:
        min = min.evalf()
        if min.is_real:
            xpoint = float(min)
            ypoint = function_numeric(xpoint)
            plt.text(xpoint, ypoint, "\u25CF", ha='center', va='center', fontsize=10)
            # plt.annotate(f'{xpoint:0.2f}', xy=(xpoint, ypoint), xytext=(xpoint - 0.1, ypoint + 0.2), fontsize=13)


    # Finish plot setup
    plt.scatter(domain, range, 0.5, marker=".")
    plt.ylim(bottom=y_bottom, top = y_top)
    plt.axvline(x=0, c="black")
    plt.axhline(y=0, c="black")
    fun_latex = sp.latex(function_to_eval)
    left_side = r'$\dot{x} = $'
    plt.title(f'Potential of {left_side} ${fun_latex}$')
    plt.legend()
    plt.subplots_adjust(bottom=0.4)
    plt.figtext(0.5, 0.1, f"Stable EQ's : ${sp.latex(mins)}$\n Unstable Eq's ${sp.latex(maxes)}$", ha='center', fontsize=12)



function_to_eval = 1 + x - x**3
plot_potentials(function_to_eval, -5, 5, -5, 5)
plt.show()


