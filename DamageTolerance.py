"""
simplistic program to calculate stuff like damage tolerance or something.



N_m is for Mx and T_m is for My

"""


import numpy as np
import math as m
import matplotlib.pyplot as plt
from MOI import Ixx, Izz, Ixx2, Izz2, centx, centz, d, b
import Loads as Loads




#Tension stress
def TensileStress(M,y,I):
    tensileStress = M * y / I
    return tensileStress

#K_l calculation
def K_i(M,y,I,a):
    K_i = TensileStress(M,y,I) * m.sqrt(m.pi*a)
    return K_i

def alpha(d):
    alpha = 2.5*10**(-3)/d
    return alpha

#f is a geometric factor    
def f(d):
    f = m.sqrt(1/(m.cos(m.pi**2*alpha(d)/360)))
    return f


#Final stress calculation
def Crackstress(M,y,I,d,a,r):
    crackStress=f(d)*K_i(M,y,I,a)/m.sqrt(r)
    return crackStress

#Main
r = 0.003

StressX = Crackstress(Loads.N_m,centz-b,Ixx,d,0.0025,r)
StressY = Crackstress(Loads.T_m,d-centx,Izz,d,0.0025,r)

finalStress = StressX + StressY

#plt.plot(Loads.span,finalStress)
#plt.plot(Loads.span,StressX)
#plt.plot(Loads.span,StressY)
#plt.show()
