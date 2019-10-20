import numpy as np
import math
import matplotlib.pyplot as plt 


##########################

def rk4SO(x_i, x_f, N, IC): #parameters: initial x,y- values, time step- h, function- f
    
    dt = (x_f - x_i) / float(N)
    x = np.arange(x_i, x_f + dt, dt) 
    y = np.zeros(int(N) + 1) #NOTE: Iterations must be integer-valued
    u = np.zeros(int(N) + 1)

    x[0], y[0], u[0] = IC #first elements
   

    #Use linspace for non-integer step-size. Linspace uses number of samples, N
    #np.arange() and range() uses step-size instead of number of samples, also inherently excuses end-point
    #look at mgrid and ogrid

    for i in range(1, N + 1):
        
        k1 = dt * f(x[i - 1], y[i - 1], u[i - 1])
        j1 = dt * g(u[i - 1])

        k2 = dt * f(x[i - 1] + 0.5 * dt, y[i - 1] + 0.5 * j1, u[i - 1] + 0.5 * k1)
        j2 = dt * g(u[i - 1] + (0.5 * k1))

        k3 = dt * f(x[i - 1] + 0.5 * dt, y[i - 1] + 0.5 * j2, u[i - 1] + 0.5 * k2)
        j3 = dt * g(u[i - 1] + (0.5 * k2))

        k4 = dt * f(x[i - 1], y[i - 1] + j3, u[i - 1] + k3)
        j4 = dt * g(u[i - 1] + 0.5 * k3)

        x[i] = x[i - 1] + dt
        y[i] = y[i - 1] + (j1 + 2.0 * (j2 + j3) + j4) / 6.0
        u[i] = u[i - 1] + (k1 + 2.0 * (k2 + k3) + k4) / 6.0

        np.savetxt('RK4SO2.txt', np.column_stack([x, y, u]), fmt='%.10f')

    print "dt = ", dt    
    return x[i], y[i], u[i]


def g(u):
    
    return u

def f(x, y, u): #just a function with variables x and y 

    return 2 

def rk4SOplot(RK4SO2):

    x, y, u = np.loadtxt('RK4SO2.txt', delimiter=' ', unpack=True)
   
    plt.plot(x, y, 'r:') #red dashed markers

    plt.axis([x[0], max(x)+max(x)/4, 0, max(y)+max(y)/4])
    plt.xlabel('t')
    plt.ylabel('x')
    plt.title('Numerical Sol using RK4 for a Quadratic Function')

    plt.show()
##########################

IC = 0, 0, 0

#k, m = 1, 1

#w = math.sqrt(k / m)

#print "w = ", w

x, y, u = rk4SO(0, 4, 100, IC)


rk4SOplot('RK4FO2.txt')
