import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import math
import pygame as py

"""
py.init()


screen = py.display.set_mode((800,600))

font = py.font.Font('freesansbold.ttf',32)

def line(x1,y1,x2,y2):
    py.draw.line(screen,(255,255,255),(x1,y1),(x2,y2),1)

def rect(x,y,L,r,g,b):
    py.draw.rect(screen,(r,g,b),(x,y,L,L))

def collisions(x,y,a):
    showscore = font.render(a, True, (255,255,255))
    screen.blit(showscore,(x,y))

#Wall
Wallx1 = 50
Wally1 = 50
Wallx2 = 50
Wally2 = 400

#Blocks
L = 50
m1x = Wallx2 + 150
m1y = Wally2 - L
m2x = Wallx2 + 450
m2y = Wally2 - L
speedm1 = 0
speedm2 = 0.1
m1 = 1
m2 = 100
count = 0
"""




def ydot(t,y):
    # right-hand side of differential equation
    return np.array([y[2],y[3],-2*w**2*y[0]+w**2*y[1]-c*y[2],w**2*y[0]-2*w**2*y[1]-c*y[3]])


def RKtrapsimp(t,y,h):
    k1 = ydot(t,y)
    k2 = ydot(t+h,y+h*k1)
    k3 = ydot(t+h/2,y+h/4*k1+h/4*k2)
    
    y1 = y + h/2*(k1+k2)
    y2 = y+h/6*k1+h/6*k2+2*h/3*k3
    e = y2-y1
    return y1,y2,e



def int_flame_rk12(y_0,t_0):
    '''
     Solve IVP:
       y' = y^2-y^3
       y(0) = y_0
       using an embedded RK1/2 pair
       (such as e.g. forward Euler/explicit trapezoid)
    '''
    t = t_0  # current time
    y = y_0  # current solution
    t_array = [t] # store times
    y_array = [y] # and solution values
    e_array = [np.array([0,0,0,0])] # estimated errors
    h = 0.1
    while t_array[-1] <= 100:
        while True:
            p = 2  # global error order for the lower order method
            y1, y2, e = RKtrapsimp(t,y,h)
            #print(e)
            if (math.sqrt(e[0]**2+e[1]**2))<tol:
                break
            # Compute the optimal time-step:
            # e   ~ O(h^(p+1))
            # tol ~ O(h_opt^(p+1))
            # h_opt = (tol/e)**(1/(p+1))*h;
            if np.isfinite((math.sqrt(e[0]**2+e[1]**2))):
                h = 0.8*(tol/(math.sqrt(e[0]**2+e[1]**2)))**(1/(p+1))*h
            else:
                h = h/2
        y = y2
        t = t+h
        # store the time/solution/error
        t_array.append(t)
        y_array.append(y)
        e_array.append(e)
        if (math.sqrt(e[0]**2+e[1]**2)) < 0.5*tol:
            h = 0.8*(tol/(math.sqrt(e[0]**2+e[1]**2)))**(1/(p+1))*h
            #or simply:
            #h = h*2
    
    plt.plot(t_array,[y_array[i][0] for i in range(len(y_array))])
    #plt.plot(t_array,[y_array[i][1] for i in range(len(y_array))])
   
    plt.show()


if __name__=='__main__':
    k = 2
    m = 4
    w = math.sqrt(k/m)
    c = 0.05
    y_0 = [-50,100,10,20]   # initial data
    t_0 = 0         # start time
    #T = 2/y_0       # end time
    tol = 1.0E-07   # tolerance for the error estimator
    int_flame_rk12(y_0,t_0)
    #avg = (t[-1]-t[0])/(len(t)-1)


