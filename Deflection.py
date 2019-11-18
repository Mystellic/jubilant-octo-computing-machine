from math import *
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate
import scipy.special as special


G=1*10**9
TableS=[]
Polar=[]

a= 1
b=1.5
c=1.2
d=1
z=5
taper = (1- (0.6*z)/11.445)    
beta= 10 #degrees
t= 0.001 #thin walled

#Polar moment of inertia variation

a= a*taper
b= b*taper
c=c*taper
d=d*taper
#Areas

Aa= a*t
Ab= b*t
Ac = c*t
Ad = d*t
#center x

Xa= a/2
Xb= t/2
Xc= t + (c/2)* cos(beta*(pi/180))
Xd= a-t/2

#center y
Ya= t/2
Yb= t+b/2
Yc= (c/2)*sin(beta*(pi/180)) +t/2 + d
Yd= t +d/2

#---------------------------------------------------------------------------------------------
#centroid
x= (Aa* Xa + Ab* Xb +Ac* Xc +Ad* Xd )/(Aa+Ab+Ac+Ad)
y= (Aa* Ya + Ab* Yb +Ac* Yc +Ad* Yd )/(Aa+Ab+Ac+Ad)


#---------------------------------------------------------------------------------------------
#moment of inertia x
Ix= (y-Ya)**2*Aa + 1/12 * t * b**3 + Ab* (y-Yb)**2


#Lift distribution
M= 1/z**2

def q(z):
    return M

def h(z):
    return E*Ix

def k(z):
    return 

i= sp.integrate.quad(q,0,z)
print(Ix)




