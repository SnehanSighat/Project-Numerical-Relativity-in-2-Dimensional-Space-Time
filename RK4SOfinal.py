import numpy as np
import matplotlib.pyplot as plt 
import math
import matplotlib.animation as animation

t_i = input('Please insert the initial time: ')
t_f = input('Please insert the final time: ')
N_t = input('Please insert the total number of temporal iterations: ')
w = 1

dt = (t_f - t_i) / float(N_t)
t = np.linspace(t_i, t_f, N_t + 1) 
x = np.zeros(int(N_t) + 1) #NOTE: Iterations must be integer-valued
u = np.zeros(int(N_t) + 1)
xtrue = np.zeros(int(N_t) + 1)

IC = t_i, 1, 0, 1 #first elements#

def rk4SOfinal(f1, f2, dt, N_t, IC): #parameters: initial x,y- values, time step- h

    #Use linspace for non-integer step-size. Linspace uses number of samples, N
    #np.arange() and range() uses step-size instead of number of samples, also inherently excuses end-point
    #look at mgrid and ogrid

    t[0], x[0], u[0], xtrue[0] = IC

    #f1_ = lambda t, x, u: np.asarray(f1(t, x, u))
    #f2_ = lambda t, x, u: np.asarray(f2(t, x, u))

    for i in range(0, N_t):
        
        k1 = dt * f1(t[i], x[i], u[i])
        j1 = dt * f2(t[i], x[i], u[i])

        k2 = dt * f1(t[i] + 0.5 * dt, x[i] + 0.5 * k1, u[i] + 0.5 * j1)
        j2 = dt * f2(t[i] + 0.5 * dt, x[i] + 0.5 * k1, u[i] + 0.5 * j1)

        k3 = dt * f1(t[i] + 0.5 * dt, x[i] + 0.5 * k2, u[i] + 0.5 * k2)
        j3 = dt * f2(t[i] + 0.5 * dt, x[i] + 0.5 * k1, u[i] + 0.5 * j2)

        k4 = dt * f1(t[i], x[i] + k3, u[i] + j3)
        j4 = dt * f2(t[i], x[i] + k3, u[i] + j3)

        t[i + 1] = t[i] + dt
        x[i + 1] = x[i] + (k1 + 2.0 * (k2 + k3) + k4) / 6.0
        u[i + 1] = u[i] + (j1 + 2.0 * (j2 + j3) + j4) / 6.0

        xtrue[i + 1] = math.cos(t[i])
 
    return t, x, u, xtrue

def f1(t, x, u): # u = '
    
    return u

def f2(t, x, u): # u' = F(t, x, u)
	
    return (w**2) * x * (-1) 

def rk4SOplotfinal(RK4SOfinal):

    t, x, u, xtrue = np.loadtxt('RK4SOfinal.txt', delimiter=' ', unpack=True)
   
    plt.plot(t, x, 'b', label = 'RK4 Sol') #red dashed markers
    plt.plot(t, xtrue, 'r--', label = 'Actual Sol')
    plt.axis([t[0], max(t)+max(t)/4, min(x), max(x)+max(x)/4])
    plt.xlabel('t')
    plt.ylabel('x')
    plt.title('Numerical Solution for 2nd ODE using RK4')
    plt.legend()
    plt.show()

t, x, u, xtrue = rk4SOfinal(f1, f2, dt, N_t, IC)

print 'For time-step, dt:', dt
np.savetxt('RK4SOfinal.txt', np.column_stack([t, x, u, xtrue]), fmt='%.10f')

rk4SOplotfinal('RK4SOfinal.txt')

fig, ax = plt.subplots()
line, = ax.plot(t, x, color='k')

def update(num, t, x, line):
    line.set_data(t[:num], x[:num])
    line.axes.axis([min(t), max(t), min(x), max(x)])
    return line,

ani = animation.FuncAnimation(fig, update, len(t), fargs=[t, x, line],
                              interval=0.1, blit=True)

ani.save('RK4SOfinale.mp4', fps=200)
plt.show()

