from math import *
import matplotlib.pyplot as plt
G=1*10**9
TableS=[]
Polar=[]


chord= 3.63
a= 0.45 * chord
b= 0.1011 * chord
c= 0.0985 * chord
d=0.45 * chord
taper = 1
betaD1= 0.458
betaD2= 0.784
betaR1= betaD1 * (pi/180) 
betaR2= betaD2 * (pi/180) 
t= 0.01 #thin walled

#Areas
Aa = a*t
Ab = b*t
Ac = c*t
Ad = d*t

    #center x
Xa= (a/2)*cos(betaR1)
Xb= t/2
Xc= a*cos(betaR1) -(t/2)
Xd= (d/2)*cos(betaR2)

    #center y
Ya= (a/2)*sin(betaR1)
Yb= (a)*sin(betaR1) + (b/2)
Yc= t+ (c/2)
Yd= t+c + (d/2)*sin(betaR2)

    

    #---------------------------------------------------------------------------------------------
#centroid
x= (Aa* Xa + Ab* Xb +Ac* Xc +Ad* Xd )/(Aa+Ab+Ac+Ad)
y= (Aa* Ya + Ab* Yb +Ac* Yc +Ad* Yd )/(Aa+Ab+Ac+Ad)


    

#---------------------------------------------------------------------------------------------
#moment of inertia x
Ix= (t*a**3*(sin(betaR1))**2)/12 + Aa* (y-Ya)**2 + 1/12 * t * b**3 + Ab* (y-Yb)**2 + 1/12 * t* c**3 + Ac * (y-Yc)**2 + (t*d**3 * (sin(betaR2))**2)/12 + Ad * (y-Yd)**2


    #polar moment of inertia
A= a*b
s= a+b+c+d


J=4*A**2/(s/t)
print(J)

print(Ix)
