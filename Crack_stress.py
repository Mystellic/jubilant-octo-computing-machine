"""

Calculate crack size or stress

"""

import math as m
import numpy as np
import matplotlib.pyplot as plt

def stress(a):
    stress = Kc/np.sqrt(m.pi*a)
    return stress


def cracksize(s):
    cracksize = (Kc*Kc)/(s*s)*m.pi
    return cracksize



#-----------------constants-----------------#

#Kc = 40e6      #Titanium
Kc = 27.5e6  #Aluminium
alu_failure_stress = 483e6
a = 0.025


#------------------Arrays------------------#

span = np.linspace(0, 10.8695, 600)
cracks = np.linspace(a, a, 600)
crack_stress = stress(cracks)
crack_safety = alu_failure_stress/crack_stress

plt.plot(span,crack_safety)
plt.show()
