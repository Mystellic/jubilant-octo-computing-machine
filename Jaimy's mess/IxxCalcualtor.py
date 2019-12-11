import pandas as pd
import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from MOI import Ixx, Izz, centx, centz


A = np.linspace(0, 0.021, 100)
t = np.linspace(0.001, 0.017, 100)
dens = 2810



tuples = []


for x in A:
    for y in t:


        #Weight
        Topplate = y * (1.5525 + 0.4 * 1.5525) * 10.8695 / 2
        Sideplate = y * (0.36 + 0.4 * 0.36) * 10.8695 / 2
        Stringers = x * 10.8695
        W = dens * (Topplate * 2 + Sideplate * 2 + Stringers)

        #Ixx
        Ixx = 0.045422*y + 0.005776*x

































'''






###


areas = []
thicknesses = []
weights = []
Ixxs = []
vs = []
for x in A:
    for y in t:


        dens = 2810
        Topplate = y * (1.5525 + 0.4 * 1.5525) * 10.8695 / 2
        Sideplate = y * (0.36 + 0.4 * 0.36) * 10.8695 / 2
        Stringers = x * 10.8695
        W = dens * (Topplate * 2 + Sideplate * 2 + Stringers)


        Ixx = 0.0065 * y + 0.03045*x
        v = (1/((71*10**9)  * Ixx))*2.37205929e+07
        areas.append(x)
        thicknesses.append(y)
        weights.append(W)
        Ixxs.append(Ixx)
        vs.append(v)

indices = []
for x in range(10000):
    if vs[x]> 2.75 and vs[x]<3.05:
        indices.append(x)


'''

