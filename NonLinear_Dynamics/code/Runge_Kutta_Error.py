import sympy as sp

# Defining a bunch of symbols
x = sp.Symbol("x")
f = sp.Function("f")(x)
f_prime = sp.diff(f, x, 1)
x_0 = sp.Symbol("x_0")
delta_t = sp.Symbol("\Delta t")
o_delta_t_5 = sp.Symbol("O(\Delta t^5)")
t = sp.Symbol('t')
t_0 = sp.Symbol('t_0')
x_of_t = sp.Function('x')('t')
x_t0 = x_of_t.subs(t, t_0)

# Generates the taylor series centered at center for the value target. Order is the number of terms
def generate_taylor(func: sp.Function, center, target, order):
    base_function = f.subs(x, center)
    taylor_expr = f.subs(x, center)
    for i in range(1, order):
        taylor_expr += (target - center)**i * sp.diff(base_function, center, i) * 1/sp.factorial(i)
    return taylor_expr


# Using taylor series to approximate each term in runge kutta
k1_approx = f.subs(x, x_0) * delta_t
k2_approx = (generate_taylor(f, x_0, x_0 + (sp.Rational(1,2) * k1_approx),4)) * delta_t
k3_approx = (generate_taylor(f, x_0, x_0 + (sp.Rational(1,2) * k2_approx), 4))* delta_t
k4_approx = (generate_taylor(f, x_0, x_0 + (k3_approx), 4)) * delta_t

# print("\n\n")
# print(f"k_1 = {sp.latex(k1_approx)}")
# print("")
# print(f"k_2 = {sp.latex(k2_approx)}")
# print("")
# print(f'k_3 = {sp.latex(sp.expand(k3_approx))}')
# print("")
# print(f'k_4 = {sp.latex(sp.expand(k4_approx))}')
# print("")

# Build the runge_kutta_approximation
runge_kutta_approx = x_0 + sp.Rational(1,6) * (k1_approx + (2 * k2_approx) + (2 * k3_approx) + k4_approx)
runge_kutta_approx = sp.expand(runge_kutta_approx)
# Remove all terms of degree higher than t^4

terms = runge_kutta_approx.as_ordered_terms()
rk_terms_filtered = [term for term in terms if term.as_poly(delta_t).degree() < 5]
runge_kutta_approx_filtered = 0
for term in rk_terms_filtered:
    runge_kutta_approx_filtered += term
runge_kutta_approx = runge_kutta_approx_filtered




#x(t)
x_0

#x'(t)

x_prime = f.subs(x, x_0)

#x''(t)
x_double_prime = sp.diff(f.subs(x, x_of_t)).subs(sp.diff(x_of_t, t), f.subs(x, x_0)).subs(x_of_t, x_0)
# x'''(t)
x_triple_prime = sp.diff(f.subs(x, x_of_t) * sp.diff(f, x).subs(x, x_of_t)).subs(sp.diff(x_of_t, t), f.subs(x, x_0)).subs(x_of_t, x_0)
# x''''(t)
f_prime_x_t = sp.diff(f, x).subs(x, x_of_t)
f_x_t = f.subs(x, x_of_t)
f_double_prime_x_t = sp.diff(f,x, 2).subs(x, x_of_t)

x_quad_prime = sp.diff(f_prime_x_t**2 * f_x_t + f_double_prime_x_t * f_x_t**2)
x_quad_prime = x_quad_prime.subs(sp.diff(x_of_t, t), f.subs(x, x_0))
x_quad_prime = x_quad_prime.subs(x_of_t, x_0)

# print(f"x(t) = {sp.latex(x_0)}")
# print(f"x'(t) = {sp.latex(x_prime)}")
# print(f"x''(t) = {sp.latex(x_double_prime)}")
# print(f"x'''(t) = {sp.latex(x_triple_prime)}")
# print(f"x''''(t) = {sp.latex(x_quad_prime)}")


# # Assemble Taylor series
# print("\n\n")
taylor_approx = x_0 + (delta_t * x_prime) + (delta_t**2) * x_double_prime * sp.Rational(1,2) + (delta_t**3) * x_triple_prime * sp.Rational(1,6) + (delta_t**4) * x_quad_prime * sp.Rational(1, 24) + o_delta_t_5


# # Simplify and compare
simplified_runge_kutta = sp.simplify(runge_kutta_approx)
simplified_taylor = sp.simplify(taylor_approx)
difference = sp.simplify(simplified_taylor - simplified_runge_kutta)

print("\n\nRunge-Kutta Approximation:")
print(sp.latex(sp.simplify(sp.collect(runge_kutta_approx, delta_t))))

print("\n\nTaylor Series Expansion:")
print(sp.latex(simplified_taylor))

# Printing final difference
print("\n\nDifference (Should be O(Δt^5)):")
print(f'E = {sp.latex(difference)}')

