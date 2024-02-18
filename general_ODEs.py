
def general_endpoint_collocation(lhs, rhs, a,  b, ua, ub, n):
    
    # solves a given ODE of up to second order and with the for lhs = rhs, u(a) = ua, u(b) = ub
    # where u denotes the solution of the ode
    # the polynomial collocation method is used at n points (currently these are evenly spaced but this will be updated to use chebyshev points)
    
    col_points = [ a + i*(b-a)/(n-1) for i in range(n)]
    x, u, du, ddu = sp.symbols('x u du ddu')
    coeff = sp.symbols('a0:{}'.format(n))
    p = coeff[0]
    for i in range(1, n):
        p += coeff[i]*x**i

    dp = sp.diff(p, x)
    ddp = sp.diff(dp, x)

    lhs = sp.sympify(lhs)
    rhs = sp.sympify(rhs)

    f = lhs - rhs
    f = f.subs([(ddu, ddp), (u, p), (du, dp)])

    eqs = [ua - p.subs(x, a)] + [f.subs(x, t) for t in col_points[1:-1]] + [ub - p.subs(x,b)]
    J = sp.matrices.Matrix([[sp.diff(eq, c) for c in coeff] for eq in eqs])

    xc = sp.matrices.Matrix(n, 1, [ 1 for i in range(n)])
    xp = sp.matrices.Matrix(n, 1, [ 0 for i in range(n)])
    er=1.0
    while er>1e-9:
        A = J
        b = sp.matrices.Matrix(n, 1, eqs)
        xp = xc
        for i in range(n):
            A = A.subs(coeff[i], xp[i])
            b = b.subs(coeff[i], xp[i])
        
        xc = xp - A.LUsolve(b) 
        er = 0.0
        for i in range(n):
            er += (xc[i] - xp[i])**2       

    return xc


def lfe_endpoint_collocation(lhs, rhs, a,  b, ua, ub, n):

    # linear finite element method collocation, similar to the global collocation but now
    # the functions used are linear within a small interval and zero everywhere else
    # the result returned is a list whos mth element is a tuple with the intercept and slope for the linear function on the mth interval

    col_points = [ a + i*(b-a)/(n-1) for i in range(n)]
    x, u, du, ddu = sp.symbols('x u du ddu')
    c = sp.symbols('c0:{}'.format(n))
    d = sp.symbols('d0:{}'.format(n))

    lhs = sp.sympify(lhs)
    rhs = sp.sympify(rhs)
    f = lhs - rhs

    # define the collocation point equations
    col_eqs = [f.subs([(u, c[i] + d[i]*col_points[i]), (du, d[i]), (ddu, 0)]) for i in range(n-1)]

    # define the continuity conditions
    cont_eqs = [c[0] + d[0]*col_points[0] - ua] \
                + [(c[i] - c[i+1]) + (d[i] - d[i+1])*col_points[i+1] for i in range(n-1)] \
                + [c[-1] + d[-1]*col_points[-1] - ub]

    # solve the system via Newtons method
    eqs = col_eqs + cont_eqs
    coeff = c + d
    J = sp.matrices.Matrix([[sp.diff(eq, c) for c in coeff] for eq in eqs])

    xc = sp.matrices.Matrix(n, 1, [ 1 for i in range(n)])
    xp = sp.matrices.Matrix(n, 1, [ 0 for i in range(n)])
    er=1.0
    while er>1e-9:
        A = J
        b = sp.matrices.Matrix(n, 1, eqs)
        xp = xc
        for i in range(n):
            A = A.subs(coeff[i], xp[i])
            b = b.subs(coeff[i], xp[i])
        
        xc = xp - A.LUsolve(b) 
        er = 0.0
        for i in range(n):
            er += (xc[i] - xp[i])**2       

    return 

a = general_endpoint_collocation('ddu', '6*x', 0, 1, 0, 1, 7)
print(a)

for m in range(3, 9, 2):
    b = general_endpoint_collocation('du', 'u', 0, 1, 1, np.e, m)
    print(b)
