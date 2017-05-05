# -*- coding: utf-8 -*-
"""
Created on Mon May  1 17:57:31 2017

@author: charlesgulian
"""

import numpy as np
import matplotlib.pyplot as plt

# Brownian Motion presentation/graphics

# Random walk generator
def random_walk_1D(T,x0=0):
    t = np.linspace(0.0,T,T+1) # Time vector
    x = np.zeros(T+1) # Position vector
    x0 = 0.0 # Initial x-value
    x[0] = x0
    steps = 2*np.random.randint(0,2,T) - 1 # Create random vector of steps of size +/- 1
    for i in range(1,len(t)):
        x[i] = x[i-1] + steps[i-1]
    return t,x,steps
    
plot_Single = False
#plot_Single = True
if plot_Single:
    time,x,steps = random_walk_1D(1000)
    plt.plot(time,x,'b')
    plt.title('Random Walk')
    plt.xlabel('Time')
    plt.ylabel('X')
    plt.show()
    
plot_Multiple = False
#plot_Multiple = True
if plot_Multiple:
    time,x1,s = random_walk_1D(1000)
    time,x2,s = random_walk_1D(1000)
    time,x3,s = random_walk_1D(1000)
    time,x4,s = random_walk_1D(1000)
    time,x5,s = random_walk_1D(1000)
    time,x6,s = random_walk_1D(1000)
    plt.plot(time,x1,'b')
    plt.plot(time,x2,'g')
    plt.plot(time,x3,'k')
    plt.plot(time,x4,'m')
    plt.plot(time,x5,'c')
    plt.plot(time,x6,'y')
    plt.title('Random Walk')
    plt.xlabel('Time')
    plt.ylabel('X')
    plt.show()
    
plot_2D = False
#plot_2D = True
if plot_2D:
    time,x,s = random_walk_1D(1000)
    time,y,s = random_walk_1D(1000)
    plt.plot(x,y,'b')
    plt.plot(0.0,0.0,'go')
    plt.plot(x[-1],y[-1],'ro')
    total_min = np.min([np.min(x),np.min(y)])
    total_max = np.max([np.max(x),np.max(y)])
    plt.axis([total_min-1,total_max+1,total_min-1,total_max+1])
    plt.grid()
    plt.title('Random Walk Path in 2D')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

# Now, let's do statistics on the random walks at different time steps

statistics = False
statistics = True
if statistics:
    X = np.zeros([10000,1001])
    for i in range(10000):
        X[i,:] = random_walk_1D(1000)[1]
    time = np.linspace(0.0,1000,1001)
    
    std_vector = np.zeros(1001)
    for i in range(1001):
        std_vector[i] = np.std(X[:,i])

    for i in range(10000):
        plt.plot(time,X[i,:],'0.8',alpha=0.10)
    plt.plot(time,std_vector,'r')
    plt.plot(time,-std_vector,'r')
        
    plt.title('Random Walk Density for 1000 Particles')
    plt.xlabel('Time')
    plt.ylabel('X')
    


