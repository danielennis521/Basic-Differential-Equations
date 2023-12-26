# some simple tools that I find helpful for solving ODE's

import numpy as np
import matplotlib.pyplot as plt


def explicit_euler(f, h, t, y0):
    return y0 + h*f(t, y0)


# implicit euler method that uses the secant method to solve the resulting nonlinear equation
def backwards_euler_secant(f, h, t, y0):
    # y0 = y1 + h*f(t, y0) => 0 = y1 - y0 + h*f(t, y1)
    # initial guesses for solution
    yp = explicit_euler(f, h, t, y0)
    yc = explicit_euler(f, h/2, t+h/2, explicit_euler(f, h/2, t, y0))
    temp = yp
    gc = (y0 - yp + h*f(t+h, yp))

    # iterate with secant method to compute solution
    while np.abs(yc - yp) > 1e-6:
        gp = gc
        gc = (y0 - yc + h*f(t+h, yc))
        temp = yc
        yc = yc - gc*(yc - yp)/(gc - gp)
        yp = temp

    return yc


# uses the second order implicit tradezoid method and then secant method to solve
def trapezoid_secant(f, h, t, y0):
    # y0 = y1 + h*(f(t, y0) + f(t, y1))/2
    # initial guesses for solution
    yp = explicit_euler(f, h, t, y0)
    yc = explicit_euler(f, h/2, t+h/2, explicit_euler(f, h/2, t, y0))
    temp = yp
    gc = (y0 - yp + h*(f(t, y0) + f(t+h, yp)))/2.0

    # iterate with secant method to compute solution
    while np.abs(yc - yp) > 1e-6:
        gp = gc
        gc = y0 - yc + h*(f(t, y0) + f(t+h, yc))/2.0
        temp = yc
        yc = yc - gc*(yc - yp)/(gc - gp)
        yp = temp

    return yc


# finds an explicit multistep scheme using n previous points of the solutions derivative
def generate_ex_multistep_scheme(n):
    # NB - the coefficients are returned in the order: a, b1, b2, b3, ... 
    #       corresponding to y_k+1 = a*y_k + b1*y'_k-n + b2*y'_k-n+1 + ...

    # construct the rhs matrix
    a = []
    for i in range(0, n):
        a.append([0, 1])

        for j in range(2, n+1):
            a[i].append(a[i][j-1] * i)
        
        for j in range(2, n+1):
            a[i][j] = a[i][j] * j

    t = [1]
    for i in range(1, n+1):
        t.append(t[i-1] * (n-1))
    a.insert(0,t)

    # construct the lhs vector
    v = [1]
    for i in range(1, n+1):
        v.append(v[i-1] * n)
    
    # solve
    print(a)
    print(v)
    a = np.array(a)
    v = np.array(v)
    print(a.transpose())
    return np.linalg.solve(a.transpose(), v)
