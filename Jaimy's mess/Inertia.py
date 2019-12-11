
import numpy as np

from math import *
import matplotlib.pyplot as plt

#t = 0.003
#AS = 0.010792
#AS_lower = 0

b = np.linspace(0.349, 0.1396, 600)
a = np.linspace(1.55, 0.620, 600)

#centroid
def Inertia_Calc(t_flange, t_spar, AS, AS_lower):
    
    centy = (a*t_flange*b/2+a*t_flange*-b/2+b*t_spar*0+b*t_spar*0+AS*b/2+AS_lower*-b/2)/(2*a*t_flange+2*b*t_spar+AS+AS_lower)
    centx = (a*t_flange*0+a*t_flange*0+b*t_spar*a/2+b*t_spar*-a/2)/(2*a*t_flange+2*b*t_spar+AS+AS_lower)
    Ixx = 2*(1/12*t_spar*b**3+b*t_spar*(centy)**2)+a*t_flange*(b/2-centy)**2+a*t_flange*(b/2+centy)**2 + AS*(b/2-centy)**2+AS_lower*(b/2+centy)**2

    Iyy = 2*(1/12*t_flange*a**3+a*t_flange*(centx)**2)+b*t_spar*(a/2-centx)**2+b*t_spar*(a/2+centx)**2+AS/2*(a/4-centx)**2+AS/2*(a/4+centx)**2+AS_lower/2*(a/4-centx)**2+AS_lower/2*(a/4+centx)**2
    return centx, centy, Ixx, Iyy

