import numpy as np
import math as m

# finds the centered finite difference scheme for the first derivative 
# the "order" parameter determines the order of accuracy
# the returned results are the coefficients a1*f(x+h) + a2*f(x-h) + a3*f(x+2h) + a4*f(x-2h) + ... = f'(x)
# in the order [a1, a2, a3, a4, ...]
# note that this could be done more simply using matricies but the method here takes advantage of the structure 
# of the problem to achieve a better run time, the number of rows we do operations on during each loop is cut in half 
def centered_first_deriv(order):
    if order%2:
        order += 1

    eqs = []
    for i in range(1, 2**(order//2 - 1) + 1):                # build differences of taylor expansions that we start with
        eqs.append([[1, -1], [2*i**(2*j+1)/m.factorial(2*j+1) for j in range(order//2)]])

    for i in range(order//2-1):

        for j in range(0, 2**(order//2 - 1), 2**(i+1)):   # weird indexing is so that we can do everything "in place"
            c = -eqs[j][1][i+1]/eqs[j+2**i][1][i+1]          # combine adjacent equations to cancel the next highest order term in each pair
            for k in range(order//2):
                eqs[j][1][k] += c*eqs[j+2**i][1][k]
            eqs[j][0] = eqs[j][0] + [c*f for f in eqs[j+2**i][0]]

    return [x/eqs[0][1][0] for x in eqs[0][0]]

print(centered_first_deriv(4))
