import numpy as np
import matplotlib.pyplot as plt 
import math
import os
import time

t_i = 0 
t_f = 5
N_t = 10000

x_i = 0
L = 1
N_x = 20

dt = (t_f - t_i) / float(N_t)
dx = (L - x_i) / float(N_x)

t = np.linspace(t_i, t_f, N_t + 1)
x = np.linspace(x_i, L, N_x + 1)
u = np.zeros((int(N_t) + 1, int(N_x) + 1))

t[0] = t_i
x[0] = x_i
u[0] = 0

def pexact(x):
	return (np.sin(2*math.pi*x)/L)

IC = np.zeros(N_x + 1)
IC[0] = 0
IC[1:] = pexact(x[1:])  #spatial initial conditions

def pde(f1, f2, dt, IC): 

    p = np.zeros((N_t + 1, N_x + 1)) #NOTICE DOUBLE PARENTHESIS
    p[0] = IC

    #f1_ = lambda t, p, u: np.asarray(f1(t, x, u))
    #f2_ = lambda t, p, u: np.asarray(f2(t, x, u))

    for i in range(0, N_t):
        
        k1 = dt * f1(t[i], p[i], u[i])
        j1 = dt * f2(t[i], p[i], u[i])

        k2 = dt * f1(t[i] + 0.5 * dt, p[i] + 0.5 * k1, u[i] + 0.5 * j1)
        j2 = dt * f2(t[i] + 0.5 * dt, p[i] + 0.5 * k1, u[i] + 0.5 * j1)

        k3 = dt * f1(t[i] + 0.5 * dt, p[i] + 0.5 * k2, u[i] + 0.5 * k2)
        j3 = dt * f2(t[i] + 0.5 * dt, p[i] + 0.5 * k1, u[i] + 0.5 * j2)

        k4 = dt * f1(t[i], p[i] + k3, u[i] + j3)
        j4 = dt * f2(t[i], p[i] + k3, u[i] + j3)

        t[i + 1] = t[i] + dt
        p[i + 1] = p[i] + (k1 + 2.0 * (k2 + k3) + k4) / 6.0
        u[i + 1] = u[i] + (j1 + 2.0 * (j2 + j3) + j4) / 6.0

    return t, p, u

def f1(t, p, u): 
    
    return u

def f2(t, p, u):

	f2 = np.zeros(N_x + 1)
	# p = np.zeros(N_x + 1)

	f2[0] = 0

	for j in range(1, N_x):

		f2[j] = (p[j + 1] - 2*p[j] + p[j-1])*(1 / dx**2)

	f2[N_x] = 0

	return f2

t, p, u = pde(f1, f2, dt, IC)	

print 'dt =', dt, 'dx =', dx

t0 = time.clock()
t1 = time.clock()
plt.ion()
y = p[0,:]
lines = plt.plot(x, y)
plt.axis([x_i, L, -2, 2])
plt.xlabel('x')
plt.ylabel('p(x,t)')
counter = 0
speed = 1
for i in range(0, p.shape[0]):
	print t[i]
	plot = True if i <= speed else i % 10 == 0
	lines[0].set_ydata(p[i,:])
	plt.legend(['t = %.0f' % t[i]])
	plt.draw()
	if plot:
		plt.savefig('kg_%04d.png' % counter)
		counter += 1
#time.sleep(0.2)