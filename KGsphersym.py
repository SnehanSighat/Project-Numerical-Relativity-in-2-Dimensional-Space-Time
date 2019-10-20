import numpy as np
import matplotlib.pyplot as plt 
import math
import time

t_i = 0 
t_f = 7
N_t = 500

r_i = 0
R = 15
N_r = 200

dt = (t_f - t_i) / float(N_t)
dr = (R - r_i) / float(N_r)

t = np.linspace(t_i, t_f, N_t + 1)
r = np.linspace(r_i, R, N_r + 1)
u = np.zeros((int(N_t) + 1, int(N_r) + 1))

t[0] = t_i
r[0] = r_i
u[0] = 0 #d_dot(0, r) = 0

ev = 5 #expectation value
var = 0.5 #variance
a = (1 / ((var*2*math.pi)**0.5)) #normalisation constant

def pprofile(r): #p(0, r)
	return a*np.exp(-((((r - ev) / ((2**0.5)*var)))**2))

IC = np.zeros(N_r + 1)
IC[0] = 0 #p(0,0)
IC[1:] = pprofile(r[1:])  #spatial initial conditions

def pde2(f1, f2, dt, IC): 

    p = np.zeros((N_t + 1, N_r + 1)) #NOTICE DOUBLE PARENTHESIS
    p[0] = IC

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

	f2 = np.zeros(N_r + 1)

	f2[0] = (2 / dr**2)*(p[1] - p[0]) #p(t, r[0] = 0)

	for j in range(1, N_r): #p(t, r[j]) excludes N_r

		f2[j] = (p[j + 1] - 2*p[j] + p[j-1])*(1 / dr**2) + (1 / (dr*r[j]))*(p[j+1] - p[j-1])

	f2[N_r] = 0 #p(t, r[N_r])

	return f2

t, p, u = pde2(f1, f2, dt, IC)	

print 'dt =', dt, 'dx =', dr

t0 = time.clock()
t1 = time.clock()
plt.ion()
y = p[0,:]
lines = plt.plot(r, y)
plt.axis([-1, R + 1, min(y) - 1, max(y) + 1])
plt.xlabel('r')
plt.ylabel('p(r,t)')
counter = 0
for i in range(0, p.shape[0]):
	print t[i]
	lines[0].set_ydata(p[i,:])
	plt.legend(['t = %.0f' % t[i]])
	plt.draw()
	plt.savefig('sphsym/kg_%04d.png' % counter)
	counter += 1