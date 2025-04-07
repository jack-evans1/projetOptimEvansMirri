import numpy as np
from scipy import optimize
from casadi import *
import time


K = 5
N = 50
L = 10
ds = L/N
alpha = 1

x0 = [ 0 for i in range(N+1)]
x0[0] = 0
x0[1] = ds
for i in range(2, N-1):
    x0[i] = i/(N-3) * (L**2 -1 -2*ds)
x0[N-1] = (L-2*ds)**2 -1 + ds
x0[N] = (L-2*ds)**2 -1 + 2*ds

y0 = [ 0 for i in range(N+1)]
y0[0] = 1
y0[1] = 1
for i in range(2, N-1):
    y0[i] = -i/(N-3)
y0[N-1] = 0
y0[N] = 0


def f(x, y, w):
    s = 0
    for i in range(N+1):
        s+= w[i]
    return K*(x[N]-L/2)**2 + ds*s


def c1e(x) :
    return -x[0]
def c2e(y) :
    return -(y[0] - 1)
def c3e(y):
    return -(y[0] - y[1])
def c4e(y):
    return -(y[N] - y[N-1])
def c5e(y):
    return -y[N]

def c1i(w, y):
    for i in range(N+1):
        return -(-w[i] + y[i])
def c2i(w, y):
    for i in range(N+1):
        return (y[i] + w[i])
def c3i(x, y):
    for i in range(N+1):
        return (1/ds*((x[i+2]-x[i+1])*(x[i+1]-x[i])+(y[i+2] - y[i+1])*(y[i+1]-y[i]))) + cos(alpha * ds)

xs_c = optimize.minimize(f, (x0, y0), method='SLSQP',
                         constraints= [{"type": "ineq", "fun":c1i}, {"type": "ineq", "fun":c2i}, {"type": "ineq", "fun":c3i}, {"type": "eq", "fun":c1e}, {"type": "eq", "fun":c2e}, {"type": "eq", "fun":c3e}, {"type": "eq", "fun":c4e}, {"type": "eq", "fun":c5e}])
print(xs_c)


# opti = Opti()
# x = opti.variable(n)
# opti.minimize(f(x))
# opti.subject_to(c1e(x))
# opti.subject_to(c2e(y))
# opti.subject_to(c3e(y))
# opti.subject_to(c4e(y))
# opti.subject_to(c1i(w, y))
# opti.subject_to(c2i(w, y))
# opti.subject_to(c3i(x, y))
# opti.solver('ipopt')
# sol = opti.solve() 