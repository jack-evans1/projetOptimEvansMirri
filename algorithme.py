import numpy as np
from scipy import optimize
from casadi import *
import time
import matplotlib.pyplot as plt


K = 5
N = 50
L = 10
ds = L/N
alpha = 1


x0 = np.zeros(N+1)
y0 = np.zeros(N+1)
theta0 = np.zeros(N)
y0[0] = 1
y0[-1] = 0
theta0[0] = 0
theta0[N-1] = 0
x0[0] = 0
x0[-1] = np.sqrt((L-2*ds)**2 -1) + 2*ds

x0[2:N-1] = np.linspace(ds, np.sqrt(L**2 -1) + ds, N-3)
y0[2:N-1] = np.linspace(1, 0, N-3)
for i in range(1, N-1):
    theta0[i] = 1/(L-2*ds)


# x0 = np.linspace(0, np.sqrt(L**2 -1), N+1)
# y0 = np.linspace(1, 0, N+1)

z0 = np.concatenate([x0, y0, np.abs(y0), theta0])



def f(z):
    s = 0
    x = z[:N+1]
    w = z[2*(N+1):]
    for i in range(N+1):
        s+= w[i]
    return K*(x[N]-L/2)**2 + ds*s


def c1e(z) :
    x = z[:N+1]
    y = z[N+1:2*(N+1)]
    w = z[2*(N+1):]
    return -x[0]
def c2e(z) :
    x = z[:N+1]
    y = z[N+1:2*(N+1)]
    w = z[2*(N+1):]
    return -(y[0] - 1)
def c3e(z):
    x = z[:N+1]
    y = z[N+1:2*(N+1)]
    w = z[2*(N+1):]
    return -(y[0] - y[1])
def c4e(z):
    x = z[:N+1]
    y = z[N+1:2*(N+1)]
    w = z[2*(N+1):]
    return -(y[N] - y[N-1])
def c5e(z):
    x = z[:N+1]
    y = z[N+1:2*(N+1)]
    w = z[2*(N+1):]
    return -y[N]
def c6e(z):
    x = z[:N+1]
    y = z[N+1:2*(N+1)]
    w = z[2*(N+1):]
    longueur = []
    for i in range(N-1):
        long = np.sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2) - ds
        longueur.append(long)
    return longueur

def c1i(z):
    x = z[:N+1]
    y = z[N+1:2*(N+1)]
    w = z[2*(N+1):]
    return -(-w + y)
def c2i(z):
    x = z[:N+1]
    y = z[N+1:2*(N+1)]
    w = z[2*(N+1):]
    return (y + w)
def c3i(z):
    x = z[:N+1]
    y = z[N+1:2*(N+1)]
    w = z[2*(N+1):]
    constraints = []
    for i in range(N - 1):
        dx1 = x[i+1] - x[i]
        dy1 = y[i+1] - y[i]
        dx2 = x[i+2] - x[i+1]
        dy2 = y[i+2] - y[i+1]
        dot = dx1 * dx2 + dy1 * dy2
        norm1 = np.sqrt(dx1**2 + dy1**2)
        norm2 = np.sqrt(dx2**2 + dy2**2)
        cos_angle = dot / (norm1 * norm2 + 1e-8)  # évite division par 0
        constraints.append(cos_angle - np.cos(alpha * ds))
    return np.array(constraints)

z = optimize.minimize(f, z0, method='SLSQP',
                         constraints= [{"type": "ineq", "fun":c1i}, {"type": "ineq", "fun":c2i}, {"type": "ineq", "fun":c3i}, {"type": "eq", "fun":c1e}, {"type": "eq", "fun":c2e}, {"type": "eq", "fun":c3e}, {"type": "eq", "fun":c4e}, {"type": "eq", "fun":c5e}, {"type": "eq", "fun":c6e}])
#print(z)

z_opt = z.x
plt.plot(z0[:N+1], z0[N+1:2*(N+1)], label='courbe initiale')
plt.plot(z_opt[:N+1], z_opt[N+1:2*(N+1)], label='courbe optimisée')
plt.legend()
plt.grid(True)
plt.title("Optimisation du tuyau de jardin")
plt.show()