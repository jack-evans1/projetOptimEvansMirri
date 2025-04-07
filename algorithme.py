import numpy as np
from scipy import optimize
from casadi import *
import time


K = 5
N = 
L = 
ds = L/N

def f(x, y, w):
    s = 0
    for i in range(N+1):
        s+= w[i]
    return K*(x[N]-L/2)**2 + ds*s


def c1e(x) :
    return 

opti = Opti()
x = opti.variable(n)
opti.minimize(f(x))
opti.subject_to(c1(x))
opti.subject_to(c2(x))
opti.solver('ipopt')
sol = opti.solve()