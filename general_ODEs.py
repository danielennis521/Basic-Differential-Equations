

# solves a given ODE of up to second order and with the for lhs = rhs, u(a) = ua, u(b) = ub
# where u denotes the solution of the ode
# the polynomial collocation method is used at n points (currently these are evenly spaced but this will be updated to use chebyshev points)
# currently Newtons method is used to solve the resulting system of equations no matter what, eventually Id like to upgrade this to have 
# a way to determine if the system is linear and change the solution appropriately
def general_endpoint_collocation(lhs, rhs, a,  b, ua, ub, n):

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
    while er>1e-6:
        A = J
        b = sp.matrices.Matrix(n, 1, eqs)
        xp = xc
        for i in range(n):
            A = A.subs(coeff[i], xp[i])
            b = b.subs(coeff[i], xp[i])
        
        xc = xp - A.LUsolve(b) 
        print(xc)

        er = 0.0
        for i in range(n):
            er += (xc[i] - xp[i])**2 
            
    return


general_endpoint_collocation('ddu', '6*x', 0, 1, 0, 1, 3)
