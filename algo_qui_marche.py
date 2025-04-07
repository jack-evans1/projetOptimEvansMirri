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
theta0[-1] = 0
x0[0] = 0
x0[-1] = np.sqrt((L-2*ds)**2 -1) + 2*ds

x0[1:N-2] = np.linspace(ds, np.sqrt(L**2 -1) + ds, N-3)
y0[1:N-2] = np.linspace(1, 0, N-3)
for i in range(1, N-1):
    theta0[i] = 1/(L-2*ds)


# x0 = np.linspace(0, np.sqrt(L**2 -1), N+1)
# y0 = np.linspace(1, 0, N+1)

z0 = np.concatenate([x0, y0, np.abs(y0), theta0])

def tuyau(theta):
    x = np.zeros(N+1)
    y = np.zeros(N+1)
    x[0] = 0
    y[0] = 1
    for i in range(N):
        x[i+1] = x[i] + ds * np.cos(theta[i])
        y[i+1] = y[i] + ds * np.sin(theta[i])
    return x, y


def f(z):
    s = 0
    theta = z[3*(N+1):]
    x, y = tuyau(theta)
    w = z[2*(N+1):3*(N+1)]
    for i in range(N+1):
        s+= w[i]
    return K*(x[N]-L/2)**2 + ds*s


def c1e(z) :
    theta = z[3*(N+1):]
    x, _ = tuyau(theta)
    return -x[0]
def c2e(z) :
    theta = z[3*(N+1):]
    _, y = tuyau(theta)
    return -(y[0] - 1)
def c3e(z):
    theta = z[3*(N+1):]
    return theta[0]
def c4e(z):
    theta = z[3*(N+1):]
    return theta[N-1]
def c5e(z):
    theta = z[3*(N+1):]
    _, y = tuyau(theta)
    return y[N] 


def ci_theta(z): #contrainte sur les angles theta
    theta = z[3*(N+1) :]
    contrainte = []
    for i in range(N-1):
        cond = alpha*ds - (theta[i+1] - theta[i])
        cond2 = alpha*ds + (theta[i+1] - theta[i])
        contrainte.append(cond)
        contrainte.append(cond2)
    return contrainte
def c1i(z):
    theta = z[3*(N+1):]
    y = tuyau(theta)[1]
    w = z[2*(N+1):3*(N+1)]
    return -(-w + y)
def c2i(z):
    theta = z[3*(N+1):]
    y = tuyau(theta)[1]
    w = z[2*(N+1):3*(N+1)]
    return (y + w)

print("z0 shape:", z0.shape)
print("should be:", 3*(N+1) + N)  # doit être 203 pour N=50
z = optimize.minimize(f, z0, method='SLSQP',
                         constraints= [{"type": "ineq", "fun":c1i}, {"type": "ineq", "fun":c2i}, {"type": "ineq", "fun":ci_theta}, {"type": "eq", "fun":c1e}, {"type": "eq", "fun":c2e}, {"type": "eq", "fun":c3e}, {"type": "eq", "fun":c4e}, {"type": "eq", "fun":c5e}])
#print(z)

z_opt = z.x
theta_opt = z_opt[3*(N+1):]
x_opt, y_opt = tuyau(theta_opt)

plt.plot(x0, y0, label='courbe initiale')
plt.plot(x_opt, y_opt, label='courbe optimisée')
plt.legend()
plt.grid(True)
plt.title("Optimisation du tuyau de jardin")
plt.show()