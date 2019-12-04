"""

Calculate crack size or stress

"""

import math as m


Kc = 27.5e6

def stress(a):
    stress = Kc/m.sqrt(m.pi*a)
    return stress


def cracksize(s):
    cracksize = (Kc*Kc)/(s*s)*m.pi
    return cracksize
