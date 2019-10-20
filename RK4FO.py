import numpy as np
import math
import matplotlib.pyplot as plt 


##########################

def rk4FO(x_i, x_f, N, IC): #parameters: initial x,y- values, time step- h, function- f
    
    h = (x_f - x_i) / float(N)
    x = np.linspace(x_i, x_f, N + 1) 
    y = np.zeros(int(N) + 1) #NOTE: Iterations must be integer-valued
    
    x[0], y[0] = IC #first elements
   

    #Use linspace for non-integer step-size. Linspace uses number of samples, N
    #np.arange() and range() uses step-size instead of number of samples, also inherently excuses end-point
    #look at mgrid and ogrid

    for i in range(1, N + 1):
        
        k1 = h * f(x[i - 1], y[i - 1])
        k2 = h * f(x[i - 1] + 0.5 * h, y[i - 1] + 0.5 * k1)
        k3 = h * f(x[i - 1] + 0.5 * h, y[i - 1] + 0.5 * k2)
        k4 = h * f(x[i - 1] + h, y[i - 1] + k3)

        x[i] = x[i - 1] + h
        y[i] = y[i - 1] + (k1 + 2.0 * (k2 + k3) + k4) / 6.0

        np.savetxt('RK4FO.txt', np.column_stack([x, y]), fmt='%.10f')

    return x[i], y[i]

def f(x, y): #just a function with variables x and y 
    
    return math.cos(x) 

def rk4FOplot(RK4FO):

    x, y = np.loadtxt('RK4FO.txt', delimiter=' ', unpack=True)
   
    plt.plot(x, y, color='black', linestyle='dashed', marker='o', markerfacecolor='red', markersize=5) #red dashed markers

    plt.axis([x[0], max(x)+max(x)/4, min(y), max(y)+max(y)/4])
    plt.xlabel('t')
    plt.ylabel('x')
    plt.title('Numerical Solution for 1st ODE using RK4')

    plt.show()
##########################

IC = (0, 0)

x, y = rk4FO(0, 20, 1000, IC)

rk4FOplot('RK4FO.txt')

